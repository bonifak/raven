"""
    xtvReader.py is a native Python module for retrieving values from a 
    TRACE XTV file.  It exposes one class that contains several methods
    which make this possible.

    This module is a bit confusing upon review, but this is mostly
    because it's built around a somewhat clunky specification for the file format.
    Basically, when an XTVFILE object is instantiated, it reads the header
    information and spawns additional objects for each component.

    The header information is sufficient for the code to determine exactly which
    bytes must be read to extract a specific data point at a specific edit.  Then
    there are a series of methods to grab either a single data point, or a list of
    points (or point pairs) along an x or z axis, or values over time.  If times
    or axial locations are requested that are between edits or mesh indices, the
    value is interpolated linearly from adjacent edits or mesh points.

    How is this module different from AptPlot and PyPost?  Well, those tools are 
    standalone applications that you interact with either using a GUI or by writing 
    separate batch scripts.  AptPlot contains its own custom batch language, while 
    PyPost uses an underlying Jython interpreter to allow batch scripts to be written
    using close to the full power of Python 2.7.  PyPost does, however, still carry 
    with it some klunky design choices of the AptPlot batch command syntax.  

    Conversely, this module can be imported into any existing native Python script, 
    giving you the ability to directly interact with an XTV file without needing to 
    spawn/fork an external process (which tends to be slow) or ensure that the script
    syntax conforms to the Jython variant of Python (rendering it less portable).

    Some of the limitations are that it does not contain many (or any, really) of the
    helper functions that the AptPlot and PyPost batch language contains.

"""

import xdrfile
from collections import namedtuple, OrderedDict
from itertools import chain
import traceback
import bisect
import re

_fineMeshVars = ['depthZrRxIn', 'depthZrRxOn', 'eCR50_46', 'hgap', 'hrfgi', 'hrfgo', 'hrfli', 'hrflo', 'hrfvi',
                'hrfvo', 'ihtfi', 'ihtfo', 'qchfi', 'qchfo', 'rftn', 'tchfi', 'tchfo', 'tmini', 'tmino', 'zht']

err_codes = {
             'TIME_UBOUND_ERR'     : "!! Error !! - Requested time point beyond last time edit",
             'TIME_LBOUND_ERR'     : "!! Error !! - Requested time point is before first time edit",
             'ERR_XRYT_INDEX'      : "!! Error !! - XR or YT index set > 0 for a variable that does not use them",
             'ERR_YT_INDEX'        : "!! Error !! - YT index set > 0 for a variable that does not use it",
             'ERR_AXRYT_INDEX'     : "!! Error !! - A, XR or YT index set > 0 for a variable that does not use them",
             'INVALID_CHANNEL1'    : "!! Error !! - The requested variable name, component type, and/or component ID are unknown",
             'INVALID_CHANNEL2'    : "!! Error !! - While decoding channel ID, got an unknown channel name or component ID",
             'INDEX_I_LBOUND_ERR'  : "!! Error !! - Invalid mesh index value - i must be > 0",
             'INDEX_J_LBOUND_ERR'  : "!! Error !! - Invalid mesh index value - j must be > 0",
             'INDEX_K_LBOUND_ERR'  : "!! Error !! - Invalid mesh index value - k must be > 0",
             'INDEX_I_UBOUND_ERR'  : "!! Error !! - Invalid mesh index value - i is too large",
             'INDEX_J_UBOUND_ERR'  : "!! Error !! - Invalid mesh index value - j is too large",
             'INDEX_K_UBOUND_ERR'  : "!! Error !! - Invalid mesh index value - k is too large",
             'AXIAL_UBOUND_ERR'    : "!! Error !! - Invalid axial height - value extends beyond last mesh point",
             'AXIAL_LBOUND_ERR'    : "!! Error !! - Invalid axial height - value comes before first mesh point",
             'AXIAL_SCALAR_ERR'    : "!! Error !! - Axial distance requested for a scalar data channel",
            }

class _Component(object):
    def __init__(self, number, compType):
        self.number = int(number)
        self.compType = compType
        self.channels = OrderedDict()
        self.dimensions = None
        self.nJuns = 0

        self.nTempl = 0
        self.templates = []

        self.nLegs = 0
        self.sidelegs = []

        self.nDynAx = 0
        self.dynAxes = []

class _Template(object):
    def __init__(self, name):
        self.name = name
        #self.index = index
        self.nCells = 0
        self.nCellsI = 0
        self.nCellsJ = 0
        self.nCellsK = 0
        self.dimi = 0
        self.dimj = 0
        self.dimk = 0
        self.coordSys = ''
        self.coordi = ''
        self.coordj = ''
        self.coordk = ''
        self.dynAxI = 0
        self.dynAxJ = 0
        self.dynAxK = 0
        self.fI = []
        self.fJ = []
        self.fK = []
        self.grav = []
        self.fa = []
        
class _SideLeg(object):
    def __init__(self, sCell, eCell, jCell):
        self.sCell = sCell
        self.eCell = eCell
        self.jCell = jCell
        
class _DynamicAxis(object):
    def __init__(self):
        self.dsAx = ''
        self.varType = ''
        self.sVarName = ''
        self.lVarName = ''
        self.vMax = 0
        
class _Channel(object):
    def __init__(self, comp, name, startIncrement):
        self.name = name
        self.comp = comp
        self.ncells = None
        self.startIncrement = startIncrement
        
class XTVError(Exception):
    """
    Custom XTV-specific exception handler (subclassed from the generic
    built-in Exception handler class.  

    This exception is raised when a data retrieval request for the XTV file
    cannot be fulfilled because the data channel (or any of its constituent
    fields) is not correct or does not otherwise exist, or the time or the 
    z location is out of bounds for the component and variable name of 
    interest.

    """

    # Use this for any XTV-specific errors like invalid data channel names, 
    # invalid time requests, etc

    def __init__(self, errmsg):
        Exception.__init__(self, errmsg)

class XtvFile(object):
    """
    Creates a new XtvFile object.

    Takes a filehandle to an open XTV file, reads & processes its header 
    information, and returns an XtvFile object ready for further operations.
    The file handle is left open for further reading.
 
    Args:
        xtvFile (file): an open file handle to the XTV file

        verbose (logical): optional argument requesting a higher level of verbosity while
        processing the header information

    Returns:
        An XtvFile object

    Example::

       import xtvReader

       # Open the XTV fiie and save the filehandle
       xtv_file = "path/to/file.xtv"
       with open(xtv_file, 'rb') as xtvFileHandle:

           # Instantiate an XtvFile object with the open filehandle.  This will read
           # in and parse the XTV header information
           xtvObj = xtvReader.XtvFile(xtvFileHandle, verbose=True)

           .....    # do some processing on the object

    """

    #  For the status, the following values apply:
    #    0 : calculation started but no data has yet been written to the XTV file
    #    1 : the calculation is in progress.  XTV file has some data but not complete
    #    2 : the calculation is complete.  The XTV file is complete
    __StartingBlock = namedtuple('StartingBlock', 'hdrString xtvMajorV xtvMinorV revNumber '
                             + 'xtvRes nUnits nComp nSVar nDVar nSChannels nDCannels dataStart '
                             + 'dataLen nPoints status spare1 spare2 spare3 fmtString unitsSys '
                             + 'sysName osString sDate sTime title')


    def __init__(self, xtvFile, verbose=False): # file handle for opened file
        """Initializes a new XtvFile instance"""

        self.xtvFile = xtvFile
        self.components = OrderedDict()
        self.times = []
        self.verbose = verbose
        self.up = xdrfile.FileUnpacker(self.xtvFile)

        try:
            self.SB = self.__StartingBlock(self.up.unpack_string(), *tuple(chain.from_iterable((tuple(
                self.up.unpack_int() for i in xrange(17)), tuple(self.up.unpack_string() for i in
                                                                           xrange(7))))))
        except:
            print "Something went wrong unpacking the Starting Block"
            print traceback.format_exc()

        if self.SB.fmtString != "MUX":
            print "Can only parse XTV files. This file is in the format " + self.SB.fmtString
            return

        if self.verbose:
            print self.SB

        recordLength = 24
        try:
            while True:
                start = self.up.get_position()
                if self.verbose:
                    print
                    print "Block starting location: " + str(start)

                # The first three values of a block are the block type,
                # revision number, and size of the block.  The data for the block
                # follows and will depend upon what type of block it is
                blockType = self.up.unpack_string()
                if self.verbose:
                    print "Block Type = " + blockType

                revision = self.up.unpack_int()  # don't need this

                jump = self.up.unpack_int()
                if self.verbose:
                    print "Blocksize = " + str(jump)


                if blockType == "GCHd":
                    compID = self.up.unpack_int()
                    if self.verbose:
                        print "   compID = " + str(compID)

                    self.up.unpack_int()

                    compType = self.up.unpack_string()
                    if self.verbose:
                        print "   compType: " + compType

                    currentComp = _Component(compID, compType)
                    self.components[(compID,compType.strip())] = currentComp

                    currentComp.title = self.up.unpack_string()
                    currentComp.dim = self.up.unpack_int()

                    currentComp.nTempl = self.up.unpack_int()
                    currentComp.templates.append(None)  # Would be convenient to make this list start at 1.  As a kluge, make the 0 element None

                    currentComp.nJuns = self.up.unpack_int()
                    currentComp.nLegs = self.up.unpack_int()

                    nSVar = self.up.unpack_int()
                    nDVar = self.up.unpack_int()
                    nVect = self.up.unpack_int()
                    nChild = self.up.unpack_int()
                    currentComp.nDynAx = self.up.unpack_int()

                elif blockType == "GD3D":
                    currentTempl = _Template("GD3D")
                    currentTempl.nCells = self.up.unpack_int()
                    currentTempl.nCellsI = self.up.unpack_int()
                    currentTempl.nCellsJ = self.up.unpack_int()
                    currentTempl.nCellsK = self.up.unpack_int()
                    currentTempl.dimi = currentTempl.nCellsI
                    currentTempl.dimj = currentTempl.nCellsJ
                    currentTempl.dimk = currentTempl.nCellsK
                    currentTempl.dynAxI = self.up.unpack_int()
                    currentTempl.dynAxJ = self.up.unpack_int()
                    currentTempl.dynAxK = self.up.unpack_int()
                    currentTempl.coordSys = self.up.unpack_string().strip()
                    if self.verbose:
                        print "  Coordsys = " + currentTempl.coordSys
                    if currentTempl.coordSys == "CART3D":
                        currentTempl.coordi = 'x'
                        currentTempl.coordj = 'y'
                        currentTempl.coordk = 'z'
                    elif currentTempl.coordSys == "CYL3D":
                        currentTempl.coordi = 'r'
                        currentTempl.coordj = 't'
                        currentTempl.coordk = 'z'

                elif blockType == "GD3A":    # Assumption is that this block will always appear right after 'GD3D'

                    currentTempl.fI = self.up.unpack_array(self.up.unpack_double)
                    currentTempl.fJ = self.up.unpack_array(self.up.unpack_double)
                    currentTempl.fK = self.up.unpack_array(self.up.unpack_double)
                    currentTempl.grav = self.up.unpack_array(self.up.unpack_double)

                    currentComp.templates.append(currentTempl)

                elif blockType == "GD2D":
                    currentTempl = _Template("GD2D")
                    currentTempl.nCells = self.up.unpack_int()
                    currentTempl.nCellsI = self.up.unpack_int()
                    currentTempl.nCellsJ = self.up.unpack_int()
                    currentTempl.dimi = currentTempl.nCellsI
                    currentTempl.dimj = currentTempl.nCellsJ
                    currentTempl.dynAxI = self.up.unpack_int()
                    currentTempl.dynAxJ = self.up.unpack_int()
                    currentTempl.coordSys = self.up.unpack_string().strip()
                    if self.verbose:
                        print "  Coordsys = " + currentTempl.coordSys
                    if currentTempl.coordSys == "CART2D":
                        currentTempl.coordi = 'x'
                        currentTempl.coordj = 'y'
                    elif currentTempl.coordSys == "CYLRZ":
                        currentTempl.coordi = 'r'
                        currentTempl.coordj = 'z'
                    elif currentTempl.coordSys == "CYLRT":
                        currentTempl.coordi = 'r'
                        currentTempl.coordj = 't'
                    elif currentTempl.coordSys == "CYLTZ":
                        currentTempl.coordi = 't'
                        currentTempl.coordj = 'z'

                elif blockType == "GD2A":  # Assumption is that this block will always appear right after 'GD2D'

                    currentTempl.fI = self.up.unpack_array(self.up.unpack_double)
                    currentTempl.fJ = self.up.unpack_array(self.up.unpack_double)
                    currentTempl.grav = self.up.unpack_array(self.up.unpack_double)

                    currentComp.templates.append(currentTempl)

                    pass

                elif blockType == "GD1D":
                    currentTempl = _Template("GD1D")
                    currentTempl.nCells = self.up.unpack_int()
                    currentTempl.nCellsI = currentTempl.nCells
                    currentTempl.dimi = currentTempl.nCells
                    currentTempl.dynAxI = self.up.unpack_int()
                    currentTempl.coordi = 'x'

                elif blockType == "GD1A":  # Assumption is that this block will always appear right after 'GD1D'

                    currentTempl.fI = self.up.unpack_array(self.up.unpack_double)
                    currentTempl.grav = self.up.unpack_array(self.up.unpack_double)
                    currentTempl.fa = self.up.unpack_array(self.up.unpack_double)

                    currentComp.templates.append(currentTempl)
                    pass

                elif blockType == "GDLg":
                    sCell = self.up.unpack_int()
                    eCell = self.up.unpack_int()
                    jCell = self.up.unpack_int()
                    currentLeg = _SideLeg(sCell, eCell, jCell)

                    if self.verbose:
                        print "  Sideleg found for component = " + str(compID)
                        print "  Sideleg starts at cell = " + str(sCell) + " and ends at cell = " + str(eCell)
                        print "  Sideleg connected to jCell = " + str(jCell)

                    currentComp.sidelegs.append(currentLeg)

                elif blockType == "DsAx":
                    dynAxis = _DynamicAxis()
                    dynAxis.dsAx = self.up.unpack_string()
                    dynAxis.varType = self.up.unpack_string()
                    dynAxis.sVarName = self.up.unpack_string()
                    dynAxis.lVarName = self.up.unpack_string()
                    dynAxis.vMax = self.up.unpack_int()

                    currentComp.dynAxes.append(dynAxis)

                elif blockType == "VARD":
                    varName = self.up.unpack_string()
                    if self.verbose:
                        print "  Variable name = " + varName
                    currentChannel = _Channel(currentComp, varName, recordLength)
                    currentComp.channels[varName.strip()] = currentChannel
                    currentChannel.varLabel = self.up.unpack_string()
                    currentChannel.uType = self.up.unpack_string()
                    currentChannel.uLabel = self.up.unpack_string()
                    currentChannel.dimPosAt = self.up.unpack_string()
                    currentChannel.freqAt = self.up.unpack_string()
                    currentChannel.cMapAt = self.up.unpack_string()
                    currentChannel.vectAt = self.up.unpack_string()
                    currentChannel.spOptAt = self.up.unpack_string()
                    currentChannel.vectName = self.up.unpack_string()
                    currentChannel.vTmpl = self.up.unpack_int()
                    currentChannel.vLength = self.up.unpack_int()
                    if currentChannel.freqAt == "TD":
                        recordLength += currentChannel.vLength * 4

                elif blockType == "DATA":
                    break

                self.up.set_position(start+jump)
        except:
            self.opened = False

        #assert recordLength == self.SB.dataLen
        #print "Record Length = " + str(recordLength)

        # Now reach into each time edit, get its time value, and save it in an array.
        for i in range(self.SB.nPoints):
            self.up.set_position(self.SB.dataStart + i*self.SB.dataLen + 20)
            self.times.append(self.up.unpack_float())

        self.opened = True

        #print self.times

        #for comp in self.components.itervalues():
            #print comp.fI

        return


    def __parseChannelString(self, channel):
        """Given a TRACE XTV channel string, parse it into its constituent elements - channel name, component id,
           and mesh index information (radial, theta, axial).  All elements returned are strings.  Conversion to
           int is left to the caller as needed"""

        # Split the variable name from the rest of the string - we'll parse that separately

        (varName, sep, id_and_mesh) = channel.rpartition('-')  # use of rpartition gracefully handles dashes in the varName part

        # The following regex look really heinous b/c it uses named groups and unsaved groups.  The goal
        # is really to just extract the component ID and mesh indices from a string that looks something
        # like this:
        #
        #     "100A11R02T01"
        #
        # and accounting for the fact that the radial and theta parts are optional
        #
        #   (?P<id>                      : starts a named group called 'id'
        #          \d+                   : match one or more digits at the start of the string
        #   )                            : end the named group ('id')
        #   (?:                          : start an unsaved grouping
        #      A(?P<axial>               : match the letter A and start a named group called 'axial'
        #                 \d\d+          : match two or more digits
        #       )                        : end the named group ('axial')
        #      (?:                       : start an unsaved grouping
        #         R(?P<radial>           : match the letter R and start a named group called 'radial'
        #                     \d\d+      : match two or more digits
        #          )                     : end the named group ('radial')
        #         (?:                    : start another unsaved group
        #            T(?P<theta>         : match the letter T and start a named group called 'theta'
        #                       \d\d+    : match two or more digits
        #             )                  : end the named group ('theta')
        #         )?                     : end the unsaved group.  It can appear zero or one time
        #      )?                        : end the unsaved group.  It can appear zero or one time
        #   )?                           : end the unsaved group.  It can appear zero or one time

        #regex = r'(?P<id>\d+)A(?P<axial>\d\d+)(?:R(?P<radial>\d\d+)(?:T(?P<theta>\d\d+))?)?'
        regex = r'(?P<id>\d+)(?:A(?P<axial>\d\d+)(?:R(?P<radial>\d\d+)(?:T(?P<theta>\d\d+))?)?)?'
        p = re.compile(regex)

        try:
            m = p.match(id_and_mesh)
            compID = m.group('id')
            a = m.group('axial')
            r = m.group('radial')
            t = m.group('theta')
        except AttributeError:
            varName = channel
            compID = 0
            a = '0'
            r = '0'
            t = '0'

        return varName, compID, a, r, t


    def __decode_channel(self, channel):
        """Decodes an XTV data channel name into its constituent pieces.  Also identifies
           what component type the requested data channel belongs to.  Numbers denoting the
           component ID and mesh indices are returned as int, not strings"""

        # Parse the channel string into its constituent parts
        (varName, id, a, r, t) = self.__parseChannelString(channel)

        # Now figure out what the component type is for the data channel
        compType = self.__getCompType(id, varName)

        #if compType is None:
        #   raise XTVError, err_codes['INVALID_CHANNEL1']

        # Convert mesh indices and comp ID to int
        id_int = int(id)
        try:
            a_int = int(a)
        except TypeError:
            a_int = None
        try:
            r_int = int(r)
        except TypeError:
            r_int = None
        try:
            t_int = int(t)
        except TypeError:
            t_int = None

        return varName, compType, id_int, a_int, r_int, t_int


    def __getCompType(self, compID, varName):
        """For a given channel identifier, retrieve the component type that
         corresponds to that channel"""

        compType = None
        tuples = self.components.keys()
        for (id,ct) in tuples:
            if id == int(compID):
                # The same compID can pair with both 'htstr' and 'htstrc' so we also need to check
                # the variable name to make sure we return the right component type string
                if self.components.get((id,ct)).channels.get(varName):
                    compType = ct
                    break
                else:
                    continue

        return compType


    def __getDimPosAt(self, id, compType, varName):
       """Private function for getting the dimPosAt attribute for a data channel of interest
          This has been isolated into its own function so that error trapping is isolated into one spot"""

       try:
           dim = self.components.get((id,compType)).channels.get(varName.strip()).dimPosAt.strip()
       except AttributeError:
           raise XTVError, err_codes['INVALID_CHANNEL1']

       return dim


    def __transform_indices(self, id, compType, varName, xr, yt, z):
        """ Given a set of indices from an APTPlot data channel (A,R,T),
            convert them to typical i,j,k indices.  Unknown values of
            the component ID, component type, or variable name will raise
            an XTVError exception"""

        try:
            dimPosAt = self.__getDimPosAt(id, compType, varName)
        except XTVError as e:
           raise e

        if dimPosAt.startswith('3'):  # strings that start with a "3" indicate this is a 3D variable
            i = xr
            j = yt
            k = z
        elif dimPosAt.startswith('2'):
            i = xr
            j = z
            k = 0
            if yt > 0:
                raise XTVError, err_codes['ERR_YT_INDEX']
        elif dimPosAt.startswith('1'):
            i = z
            j = 0
            k = 0
            if yt > 0 or xr > 0:
                raise XTVError, err_codes['ERR_XRYT_INDEX']
        elif dimPosAt.startswith('0'):
            i = 0
            j = 0
            k = 0
            if yt > 0 or xr > 0 or z > 0:
                raise XTVError, err_codes['ERR_AXRYT_INDEX']
        else:
            print "Programming Error - __transform_indices() has encountered a variable with a dimension it can't handle"
            raise ValueError

        return i, j, k


    def __interpolate(self, x1, y1, x2, y2, x):
        return (y2 - y1) * (x-x1)/(x2-x1) + y1


    def getAxialDataChannel(self, time, channel, zLoc):
        """
        Retrieves a single XTV value for a single data channel at a particular time and
        at a particular axial (z) location
 
        Args:
           time (float): the time point from which to retrieve a value.

           channel (str): a string of the XTV data channel of interest, i.e. 'rftn-140A03R05'.

           zLoc (float): the axial location at which to retrieve the data value. This 
           implies that the 'A' field in the channel string is irrelevant'. 

        Returns:
           a single floating point value

        Raises:
           XTVError will be raised if the requested data channel or z location is not
           available in the XTV file.

        Example::

           import xtvReader

           # Open the XTV fiie and save the filehandle
           xtv_file = "path/to/file.xtv"
           with open(xtv_file, 'rb') as xtvFileHandle:

               # Instantiate an XtvFile object with the open filehandle.  This will read
               # in and parse the XTV header information
               xtvObj = xtvReader.XtvFile(xtvFileHandle, verbose=True)

               # Retrieve a single value at a particular axial location
               value = xtvObj.getAxialDataChannel(10.0, 'rftn-140A01R08', 2.45)

        """
        try:
            (varName, compType, id, a, r, t) = self.__decode_channel(channel)
            value = self.getAxialData(time, id, compType, varName, zLoc, r, t)
        except XTVError as e:
            raise e    # Don't trap any XTV-specific errors here - let the caller handle them

        return value


    def getDataChannel(self, time, channel):
        """
        Retrieves a single value at a particular time for a particular mesh index location

        Args:
           time (float): the time point from which to retrieve a value.

           channel (str): a string of the XTV data channel of interest, i.e. 'rftn-140A03R05'.

        Returns:
           a single floating point value

        Raises:
           XTVError will be raised if the requested data channel or z location is not
           available in the XTV file.


        Example::

           import xtvReader

           # Open the XTV fiie and save the filehandle
           xtv_file = "path/to/file.xtv"
           with open(xtv_file, 'rb') as xtvFileHandle:

               # Instantiate an XtvFile object with the open filehandle.  This will read
               # in and parse the XTV header information
               xtvObj = xtvReader.XtvFile(xtvFileHandle, verbose=True)

               # Retrieve a single value at a particular point in time
               value = xtvObj.getDataChannel(50.0, 'pn-55A06')

        """

        try:
            (varName, compType, id, a, r, t) = self.__decode_channel(channel)
            (i,j,k) = self.__transform_indices(id, compType, varName, r, t, a)
            value = self.getData(time, id, compType, varName, i, j, k)
        except XTVError as e:
            raise e      # Don't trap any XTV-specific errors here - let the caller handle them

        return value


    def getAxialData(self, time, id, compType, varName, zLoc, xr=0, yt=0):
        """
        The purpose of this routine is to retrieve the value for a particular data channel
        at a particular axial location at a specific point in time

        Args:
           time (float): the time point from which to retrieve a value.

           id (int): the component number from which to retrieve a value

           compType (str): the component type from which to retrieve a value

           varName (str): the XTV variable name to retrieve

           zLoc (float): the axial location at which to retrieve the data value.

           xr (int): the radial/x coordinate index of the variable given by varName.  If 
           varName does not have an X/R index, then it can be omitted.

           yt (int): the theta/y coordinate index of the variable given by varName.  If 
           varName does not have a Y/T index, then it can be omitted.

        Returns:
           a single floating point value

        Raises:
           XTVError will be raised if the requested data channel or z location is not
           available in the XTV file.

        Example::

           import xtvReader

           # Open the XTV fiie and save the filehandle
           xtv_file = "path/to/file.xtv"
           with open(xtv_file, 'rb') as xtvFileHandle:

               # Instantiate an XtvFile object with the open filehandle.  This will read
               # in and parse the XTV header information
               xtvObj = xtvReader.XtvFile(xtvFileHandle, verbose=True)

               # Retrieve a single value at a particular axial location & point in time
               value = xtvObj.getAxialData(10.0, 140, 'htstrc', 'rftn', 2.45, 8)

        """

        # Knowing the z location we want to get a value for, we first need
        # to retrieve the i indices that bound that axial location at the two
        # time points that bound the desired time.  To do that, we will use the FI array.
        # It contains the z locations of each face.
        # If the variable is an edge variable, then we just use those values to interpolate.
        # If the variable is a cell center variable, then we first need to compute the cell center z values.
        # One question - what do we do for side legs in TEEs?  Looks like iFaces is just a construction from the dx array.


        # When the user selects a time value less than zero, that is a cue to use
        # the last time point, whatever it may be
        if time < 0:
            time = self.times[-1]
        elif time > self.times[-1]:      # caller has requested a time point beyond the last edit.
            raise XTVError, err_codes['TIME_UBOUND_ERR']
        elif time < self.times[0]:      # caller has requested a time point before the first edit.
            raise XTVError, err_codes['TIME_LBOUND_ERR']

        time_index = bisect.bisect_right(self.times, float(time))
        if time == self.times[time_index-1]:

            # First, we need to get the vector of axial heights
            zht = self.getAxialLocations(time, id, compType, varName)

            # Do some input checking
            if len(zht) == 0:   # list is empty.  Must be a scalar data channel
                raise XTVError, err_codes['AXIAL_SCALAR_ERR']
            elif zLoc > zht[-1]:
                raise XTVError, err_codes['AXIAL_UBOUND_ERR']
            elif zLoc < zht[0]:
                raise XTVError, err_codes['AXIAL_LBOUND_ERR']

            # Get the lower index that bounds the requested height
            z = bisect.bisect_right(zht, float(zLoc))

            # If the requested height exactly matches axial location at z, then
            # we can use that index to retrieve the parameter of interest directly
            if zLoc == zht[z-1]:
                (i,j,k) = self.__transform_indices(id, compType,varName, xr, yt, z)
                value = self.getData(time, id, compType,varName,i,j,k)
            else:  # the requested height is between two axial levels.  Need to interpolate
                (i,j,k) = self.__transform_indices(id, compType,varName, xr, yt, z)
                lVal = self.getData(time, id, compType,varName,i,j,k)

                (i,j,k) = self.__transform_indices(id, compType,varName, xr, yt, z+1)
                uVal = self.getData(time, id, compType,varName,i,j,k)
                value = self.__interpolate(zht[z-1], lVal, zht[z], uVal, zLoc)

        else:
            timeLower = self.times[time_index-1]
            timeUpper = self.times[time_index]

            # Get the heights at the lower time bound
            zht = self.getAxialLocations(timeLower, id, compType, varName)

            # Do some input checking
            if len(zht) == 0:     # list is empty.  Must be a scalar data channel
                raise XTVError, err_codes['AXIAL_SCALAR_ERR']
            elif zLoc > zht[-1]:
                raise XTVError, err_codes['AXIAL_UBOUND_ERR']
            elif zLoc < zht[0]:
                raise XTVError, err_codes['AXIAL_LBOUND_ERR']

            # Get the lower index that bounds the requested height
            z = bisect.bisect_right(zht, float(zLoc))

            # If the requested height exactly matches the axial location at z, then
            # we can use that index to retrieve the parameter of interest directly
            if zLoc == zht[z-1]:
                (i,j,k) = self.__transform_indices(id, compType,varName, xr, yt, z)
                lvalue = self.getData(timeLower, id, compType, varName, i, j, k)
            else:  # the requested height is between two axial levels.  Need to interpolate
                (i,j,k) = self.__transform_indices(id, compType, varName, xr, yt, z)
                lVal = self.getData(timeLower, id, compType, varName, i, j, k)

                (i,j,k) = self.__transform_indices(id, compType, varName, xr, yt, z+1)
                uVal = self.getData(timeLower, id, compType, varName, i, j, k)
                lvalue = self.__interpolate(zht[z-1], lVal, zht[z], uVal, zLoc)

            # Get the heights at the upper time bound
            zht = self.getAxialLocations(timeUpper, id, compType, varName)

            # Get the lower index that bounds the requested height
            z = bisect.bisect_right(zht, float(zLoc))

            # If the requested height exactly matches the axial location at z, then
            # we can use that index to retrieve the parameter of interest directly
            if zLoc == zht[z-1]:
                (i,j,k) = self.__transform_indices(id, compType,varName, xr, yt, z)
                uvalue = self.getData(timeUpper, id, compType, varName, i, j, k)
            else:  # the requested height is between two axial levels.  Need to interpolate
                (i,j,k) = self.__transform_indices(id, compType,varName, xr, yt, z)
                lVal = self.getData(timeUpper, id, compType, varName, i, j, k)

                (i,j,k) = self.__transform_indices(id, compType,varName, xr, yt, z+1)
                uVal = self.getData(timeUpper, id, compType, varName, i ,j ,k)
                uvalue = self.__interpolate(zht[z-1], lVal, zht[z], uVal, zLoc)

            # Now perform the final interpolation at the requested time.
            value = self.__interpolate(timeLower, lvalue, timeUpper, uvalue, time)

        return value


    def getAxialLocations(self, time, id, compType, varName):
        """
        The purpose of this routine is to retrieve the axial locations that 
        correspond to a particular data channel for a given component at a 
        particular point in time.  For all but fine mesh variables (whose
        node heights can vary with time), the time value is irrelevant and not
        used.

        Args:
           time (float): the time point from which to retrieve a value.

           id (int): the component number from which to retrieve a value

           compType (str): the component type from which to retrieve a value

           varName (str): the XTV variable name associated with the axial locations that
           should be retrieved.  Edge-based variables will retrieve face locations and
           cell-centered values will retrieve cell center locations.

        Returns:
           a list of floating point values

        Raises:
           XTVError will be raised if the requested data channel or z location is not
           available in the XTV file.

        Example::

           import xtvReader

           # Open the XTV fiie and save the filehandle
           xtv_file = "path/to/file.xtv"
           with open(xtv_file, 'rb') as xtvFileHandle:

               # Instantiate an XtvFile object with the open filehandle.  This will read
               # in and parse the XTV header information
               xtvObj = xtvReader.XtvFile(xtvFileHandle, verbose=True)

               # Retrieve the fine mesh zht array for a htstr at 10.0 secs
               zht = xtvObj.getAxialLocations(10.0, 140, 'htstrc', 'rftn')

               # Retrieve the axial cell-center heights for a VESSEL cell-center variable
               cellHeights = xtvObj.getAxialLocations(0.0, 140, 'vessel', 'pn')

               # Retrieve the axial face heights for a VESSEL edge variable
               faceHeights = xtvObj.getAxialLocations(0.0, 140, 'vessel', 'vlnz')

        """

        #  There are few nuances.  For some data channels, the axial locations will correspond to
        #  cell centers, and for others, they will correspond to cell edges.

        #  For TEE-based components, the vector retrieved will include both the main tube and side tube.  It will be
        #  up to the calling routine to sort out what it wants

        #  For HTSTR components, some variables may be 2D.  For those variables, this routine will
        #  retrieve a 1D array of axial elevations/locations.  Some variables may be fine mesh variables.  In those
        #  cases, this routine will retrieve a vector of zht values that is only as large as needed.

        # For VESSEL components, only the z-direction elevation information will be retrieved.

        # If the variable is a 0D variable, then an empty list is returned.

        try:
            dimPosAt = self.__getDimPosAt(id, compType, varName)
        except XTVError as e:
           raise e

        tmplInd = self.components.get((id,compType)).channels.get(varName.strip()).vTmpl

        if dimPosAt == '0D':  # If the variable is a scalar, the entire notion of axial
            return []           # locations makes no sense, so just return an empty list

        zLocs = []
        if compType == 'htstrc':

            if dimPosAt == '1dFa':
                vLength = self.components.get((id,compType)).channels.get(varName).vLength
            elif dimPosAt == '2dFaJ':
                vLength = self.components.get((id,compType)).templates[tmplInd].dimj + 1
            else:
                print "Programming Error - unknown dimension for variable " + str(varName)
                exit()

            if varName in _fineMeshVars:
                cell = 1  # Need the entire array so start at the first mesh location
                startingPoint = self.components.get((id,compType)).channels.get('zht').startIncrement + (cell-1) * 4
                startingEdit = bisect.bisect_right(self.times, float(time))
                self.up.set_position(self.SB.dataStart + (startingEdit-1)*self.SB.dataLen + startingPoint)
                zLocs = self.up.unpack_farray(vLength, self.up.unpack_float)
                zLocs = zLocs[0:zLocs.index(-1.0)]  # Filter out all the values that are unused
            else:
                zLocs = self.components.get((id,compType)).templates[tmplInd].fJ

        elif compType == 'htstr':
            zLocs = self.components.get((id,compType)).templates[tmplInd].fI[1:]

        elif  compType == 'vessel':
            zLocs = self.components.get((id,compType)).templates[tmplInd].fK
            if dimPosAt != '3dFaK':
                zLocs = [(zLocs[i]+zLocs[i+1])*0.5 for i in range(len(zLocs)-1)]   # Transform edge locations to cell center locations by calculating midpoints

        else:
            if dimPosAt.startswith('1'):
                zLocs = self.components.get((id,compType)).templates[tmplInd].fI
                if dimPosAt == '1dCc':
                    zLocs = [(zLocs[i]+zLocs[i+1])*0.5 for i in range(len(zLocs)-1)]   # Transform edge locations to cell center locations by calculating midpoints
            else:
                print "Programming Error - Not sure how to get axial locations for the requested variable " +  str(varName)
                exit()

        zLocs = [round(z,13) for z in zLocs]   # make sure floating point precision issues don't give us weird numbers
        return zLocs


    def getData(self, time, id, compType, varName, i=1, j=0, k=0):

        """
        Retrieves a value from the XTV file for a single point in
        time for a requested component ID, component type, & variable name. It
        works by first calculating the stride & file pointer offset for the
        variable of interest.  It then uses that offset to grab the value
        directly.

        Args:
           time (float): the time point from which to retrieve a value.

           id (int): the component number of the desired data channel

           compType (str): the component type of the desired data channel

           varName (str): the XTV variable name to retrieve

           i (int): the i-coordinate index of the data channel to be retrieved.  If
           varName does not have an i index, then it can be omitted.  In that case,
           a default value of 1 will be set.

           j (int): the j-coordinate index of the data channel to be retrieved.  If
           varName does not have a j index, then it can be omitted.

           k (int): the k-coordinate index of the data channel to be retrieved.  If
           varName does not have a k index, then it can be omitted.

        Returns:
           a single floating point value

        Raises:
           XTVError will be raised if the requested data channel or z location is not
           available in the XTV file.

        Example::

           import xtvReader

           # Open the XTV fiie and save the filehandle
           xtv_file = "path/to/file.xtv"
           with open(xtv_file, 'rb') as xtvFileHandle:

               # Instantiate an XtvFile object with the open filehandle.  This will read
               # in and parse the XTV header information
               xtvObj = xtvReader.XtvFile(xtvFileHandle, verbose=True)

               # Retrieve a single value at a particular point in time
               value = xtvObj.getData(10.0, 200, 'vessel', 'vln', 2, 3, 4)

        """

        try:
            dimPosAt = self.__getDimPosAt(id, compType, varName)
        except XTVError as e:   # if the component ID, component type or variable name are unknown, raise an exception
           raise e

        # Retrieve the dimensions associated with the current variable name.  If the variable has no
        # dimensions, it is likely a '0D' varName and we won't need them.  In that case, move on like
        # nothing happened.
        tmplInd = self.components.get((id,compType)).channels.get(varName.strip()).vTmpl
        vLength = self.components.get((id,compType)).channels.get(varName.strip()).vLength
        try:
            dimi = self.components.get((id,compType)).templates[tmplInd].dimi
            dimj = self.components.get((id,compType)).templates[tmplInd].dimj
            dimk = self.components.get((id,compType)).templates[tmplInd].dimk
        except AttributeError:
            dimi = None
            dimj = None
            dimk = None

        if dimPosAt.startswith('3'):  # strings that start with a "3" indicate this is a 3D variable

            if i is None or i == 0:
                raise XTVError, err_codes['INDEX_I_LBOUND_ERR']
            if j is None or j == 0:
                raise XTVError, err_codes['INDEX_J_LBOUND_ERR']
            if k is None or k == 0:
                raise XTVError, err_codes['INDEX_K_LBOUND_ERR']

            coordSystem = self.components.get((id,compType)).templates[tmplInd].coordSys
            if dimPosAt.endswith("I"):  # strings that end with "I" correspond to XR face vectors
                nmesh_i = vLength/(dimj*dimk)
                nmesh_j = vLength/((dimi+1)*dimk)
                nmesh_k = vLength/((dimi+1)*dimj)
                levdim = (dimi+1) * dimj
                # the way the 3D cells are ordered in XTV is different from the way we typically count & label
                # them.  The calculation below reflects the fact that the radial index changes most rapidly,
                # whereas normal numbering would start counting by theta direction first, then by ring.
                cell = (k-1)*levdim + (j-1) * (dimi+1) + i
            elif dimPosAt.endswith("J"):  # strings that end with "J" correspond to YT face vectors
                if coordSystem == 'CYL3D':
                    nmesh_i = vLength/(dimj*dimk)
                    nmesh_j = vLength/(dimi*dimk)
                    nmesh_k = vLength/(dimi*dimj)
                    levdim = dimi * (dimj)
                else:
                    nmesh_i = vLength/((dimj+1)*dimk)
                    nmesh_j = vLength/(dimi*dimk)
                    nmesh_k = vLength/(dimi*(dimj+1))
                    levdim = dimi * (dimj+1)
                # the way the 3D cells are ordered in XTV is different from the way we typically count & label
                # them.  The calculation below reflects the fact that the radial index changes most rapidly,
                # whereas normal numbering would start counting by theta direction first, then by ring.
                cell = (k-1)*levdim + (j-1) * dimi + i
            elif dimPosAt.endswith("K"):  # strings that end with "K" correspond to Z face vectors
                nmesh_i = vLength/(dimj*(dimk+1))
                nmesh_j = vLength/(dimi*(dimk+1))
                nmesh_k = vLength/(dimi*dimj)
                levdim = dimi * dimj
                cell = (k-1)*levdim + (j-1) * dimi + i
            elif dimPosAt.endswith("c"):  # strings that end with "c" correspond to cell center
                nmesh_i = vLength/(dimj*dimk)
                nmesh_j = vLength/(dimi*dimk)
                nmesh_k = vLength/(dimi*dimj)
                levdim = dimi * dimj
                cell = (k-1)*levdim + (j-1) * dimi + i
            else:
                raise XTVError, "!! Programming Error !! - Unexpected dimPosAt string"

            if i > nmesh_i:
                raise XTVError, err_codes['INDEX_I_UBOUND_ERR']
            if j > nmesh_j:
                raise XTVError, err_codes['INDEX_J_UBOUND_ERR']
            if k > nmesh_k:
                raise XTVError, err_codes['INDEX_K_UBOUND_ERR']

        elif dimPosAt.startswith('2'):
            nmesh_i = vLength/(dimj+1)
            if i is None or i == 0:
                raise XTVError, err_codes['INDEX_I_LBOUND_ERR']
            elif i > nmesh_i:
                raise XTVError, err_codes['INDEX_I_UBOUND_ERR']

            nmesh_j = vLength/(dimi)
            if j is None or j == 0:
                raise XTVError, err_codes['INDEX_J_LBOUND_ERR']
            elif j > nmesh_j:
                raise XTVError, err_codes['INDEX_J_UBOUND_ERR']

            cell = (j-1) * dimi + i
        elif dimPosAt.startswith('1'):
            if i is None or i == 0:
                raise XTVError, err_codes['INDEX_I_LBOUND_ERR']
            elif i > vLength:
                raise XTVError, err_codes['INDEX_I_UBOUND_ERR']

            cell = i
        elif dimPosAt.startswith('0'):
            cell = 1
        else:
            raise XTVError, "!! Programming Error !! - Encountered a variable with a dimension I can't handle"

        startingPoint = self.components.get((id,compType)).channels.get(varName).startIncrement + (cell-1) * 4

        # When the user selects a time value less than zero, that is a cue to use
        # the last time point, whatever it may be
        if time < 0:
            time = self.times[-1]
        elif time > self.times[-1]:      # caller has requested a time point beyond the last edit.
            raise XTVError, err_codes['TIME_UBOUND_ERR']
        elif time < self.times[0]:      # caller has requested a time point before the first edit.
            raise XTVError, err_codes['TIME_LBOUND_ERR']

        # If the requested time lines up exactly with the time of one of the
        # existing graphics edits, then just grab the value directly and get out.
        startingEdit = bisect.bisect_right(self.times, float(time))
        if time == self.times[startingEdit-1]:
            self.up.set_position(self.SB.dataStart + (startingEdit-1)*self.SB.dataLen + startingPoint)
            value = self.up.unpack_float()
            return value

        # Otherwise, grab the values at time points that bound the requested time and interpolate.
        self.up.set_position(self.SB.dataStart + (startingEdit-1)*self.SB.dataLen + startingPoint)
        y1 = self.up.unpack_float()
        self.up.set_position(self.SB.dataStart + startingEdit*self.SB.dataLen + startingPoint)
        y2 = self.up.unpack_float()
        value = self.__interpolate(self.times[startingEdit-1], y1, self.times[startingEdit], y2, time)
        return value


    def getTimeData(self, times, id, compType, varName, i=1, j=0, k=0):
        """
        Given a list of time points, retrieve a list of values from the XTV 
        file for the requested component ID, component type, variable name,
        and mesh indices.

        Args:
           times (float): a list of the time points for which to retrieve 
           an XTV variable.

           id (int): the component number of the desired data channel

           compType (str): the component type of the desired data channel

           varName (str): the XTV variable name to retrieve

           i (int): the i-coordinate index of the data channel to be retrieved.  If
           varName does not have an i index, then it can be omitted.  In that case,
           a default value of 1 will be set.

           j (int): the j-coordinate index of the data channel to be retrieved.  If
           varName does not have a j index, then it can be omitted.

           k (int): the k-coordinate index of the data channel to be retrieved.  If
           varName does not have a k index, then it can be omitted.

        Returns:
           a list of floats

        Raises:
           XTVError will be raised if the requested data channel or z location is not
           available in the XTV file.

        Example::

           import xtvReader

           # Open the XTV fiie and save the filehandle
           xtv_file = "path/to/file.xtv"
           with open(xtv_file, 'rb') as xtvFileHandle:

               # Instantiate an XtvFile object with the open filehandle.  This will read
               # in and parse the XTV header information
               xtvObj = xtvReader.XtvFile(xtvFileHandle, verbose=True)

               # Retrieve a list of values for a given list of time points
               time_points = [10.0, 20.0, 30.0]
               value = xtvObj.getTimeData(time_points, 200, 'vessel', 'vln', 2, 3, 4)

        """


        result = []
        for time in times:
            try:
                result.append(self.getData(time,id,compType,varName,i,j,k))
            except XTVError as e:
                raise e

        return result


    def getTimeVector(self, channel):
        """
        Returns a list of tuples of all the (time, value) pairs for a particular data channel

        Args:
           channel (str): the XTV data channel of interest, i.e. 'rftn-140A03R05'.

        Returns:
           a list of tuples of (time, value) pairs.  Values in tuples are floats.

        Raises:
           XTVError will be raised if the requested data channel is not available in the XTV file.

        Example::

           import xtvReader

           # Open the XTV fiie and save the filehandle
           xtv_file = "path/to/file.xtv"
           with open(xtv_file, 'rb') as xtvFileHandle:

               # Instantiate an XtvFile object with the open filehandle.  This will read
               # in and parse the XTV header information
               xtvObj = xtvReader.XtvFile(xtvFileHandle, verbose=True)

               # Retrieve a list of (time, value) pairs for a particular data channel
               list = xtvObj.getTimeVector('rftn-140A01R08')

        """

        (varName, compType, id, a, r, t) = self.__decode_channel(channel)
        (i,j,k) = self.__transform_indices(id, compType, varName, r, t, a)
        vector = []
        for time in self.times:
            value = self.getData(time, id, compType, varName, i, j, k)
            vector.append((time,value))
        return vector


    def getTimeVectorAxial(self, channel, zLoc):
        """
        Returns a list of tuples of all the (time, value) pairs for a particular 
        data channel at a specific axial location.

        Args:
           channel (str): a string of the XTV data channel of interest, i.e. 'rftn-140A03R05'.

           zLoc (float): the axial location at which to retrieve the data value. This
           implies that the 'A' field in the channel string is irrelevant.

        Returns:
           a list of tuples of (time, value) pairs.  Values in tuples are floats.

        Raises:
           XTVError will be raised if the requested data channel or z location is not
           available in the XTV file.

        Example::

           import xtvReader

           # Open the XTV fiie and save the filehandle
           xtv_file = "path/to/file.xtv"
           with open(xtv_file, 'rb') as xtvFileHandle:

               # Instantiate an XtvFile object with the open filehandle.  This will read
               # in and parse the XTV header information
               xtvObj = xtvReader.XtvFile(xtvFileHandle, verbose=True)

               # Retrieve a list of (time, value) pairs at a particular axial location
               list = xtvObj.getTimeVectorAxial('rftn-140A01R08', 4.5)

        """

        (varName, compType, id, a, xr, yt) = self.__decode_channel(channel)

        vector = []
        for time in self.times:

            # First, we need to get the vector of axial heights
            zLocs = self.getAxialLocations(time, id, compType, varName)

            # Get the lower index that bounds the requested height
            z = bisect.bisect_right(zLocs, float(zLoc))

            # If the requested height exactly matches axial location at z, then
            # we can use that index to retrieve the parameter of interest directly
            if zLoc == zLocs[z-1]:
                (i,j,k) = self.__transform_indices(id, compType, varName, xr, yt, z)
                value = self.getData(time, id, compType, varName, i, j, k)
            else:  # the requested height is between two axial levels.  Need to interpolate
                (i,j,k) = self.__transform_indices(id, compType, varName, xr, yt, z)
                lVal = self.getData(time, id, compType, varName, i, j, k)

                (i,j,k) = self.__transform_indices(id, compType,varName, xr, yt, z+1)
                uVal = self.getData(time, id, compType, varName, i, j, k)
                value = self.__interpolate(zLocs[z-1], lVal, zLocs[z], uVal, zLoc)

            vector.append((time,value))

        return vector


    def getAxialVector(self, time, channel):
        """
        Returns a list of tuples of all the (axial location, value) pairs for
        a particular data channel at a particular time

        Args:
           time (float): the time point from which to retrieve a value.

           channel (str): a string of the XTV data channel of interest, i.e. 'rftn-140A03R05'.
           In this case, the "A" index has no relevance to the values retrieved although error
           processing may still require it to have a valid value.

        Returns:
           a list of tuples of (z, value) pairs.  Values in tuples are floats.

        Raises:
           XTVError will be raised if the requested data channel or z location is not
           available in the XTV file.

        Example::

           import xtvReader

           # Open the XTV fiie and save the filehandle
           xtv_file = "path/to/file.xtv"
           with open(xtv_file, 'rb') as xtvFileHandle:

               # Instantiate an XtvFile object with the open filehandle.  This will read
               # in and parse the XTV header information
               xtvObj = xtvReader.XtvFile(xtvFileHandle, verbose=True)

               # Retrieve a list of (z location, value) pairs at a particular time
               list = xtvObj.getAxialVector(10.0, 'rftn-140A01R08')
        """

        (varName, compType, id, a, r, t) = self.__decode_channel(channel)
        #(i,j,k) = self.__transform_indices(id, compType, varName, r, t, a)
        zLocs = self.getAxialLocations(time, id, compType, varName)
        vector = []
        for z in zLocs:
            value = self.getAxialData(time, id, compType, varName, z, r, t)
            vector.append((z,value))
        return vector


    #def getXData(self, time, id, compType, varName, i_s, j=0, k=0):
    #    result = []
    #    for i in i_s:
    #        result.append(self.getData(time,id,compType,varName, i, j, k))
    #    return result


    #def getZData(self, time, id, compType, varName, i=1, j=1, k_s=(0,)):
    #    result = []
    #    for k in k_s:
    #        result.append(self.getData(time,id,compType,varName, i, j, k))
    #    return result














