#!python
import numpy as np
from scipy.spatial import distance
from xml.etree.ElementTree import ElementTree as ET
from xml.etree.ElementTree import Element,SubElement,dump
from optparse import OptionParser as OP
import vtk

color_table = {
    'red':[255,0,0],
    'green':[0, 255, 0],
    'blue':[0, 0, 255],
    'white':[255, 255, 255],
    'sandybrown':[244,164,96],
    'yellow':[255,255,0],
}

class LignumVTKSpline:
    def __init__(self,branch_lss,cylinder:bool,color_name:str,radius_name:str,lgatype_name:str,resolution:int,
                 nsides:int):
        self.point_index = -1
        self.cylinder=cylinder
        #List of list of branches
        self.branch_lss=branch_lss
        self.color_name=color_name
        self.radius_name=radius_name
        self.lgatype_name=lgatype_name
        self.nsides = nsides
        self.resolution = resolution
        self.spline_points=0
        self.branch_length=0
        self.spline_segment_length=0
        #List of actors, one for each branch. Will be added to writer for file output.
        self.renderer = vtk.vtkRenderer()
        #vtkXMLMultiBlockDataWriter writes the output to ParaView compatible '.vtp' file 
        self.writer = vtk.vtkXMLMultiBlockDataWriter()
    def reset_index(self):
        self.point_index = -1
    def next_index(self):
        self.point_index = self.point_index+1
        return self.point_index
    def tree_to_tube_actors(self):
        for ts_ls in self.branch_lss:
            self.branch_to_tube_actor(ts_ls)
    def branch_to_tube_actor(self,ts_ls):
        lgm_xml=LignumXML()
        points = vtk.vtkPoints()
        #Calculate total branch length
        self.set_branch_length(ts_ls)
        i=0
        for ts in ts_ls:
            p = lgm_xml.ts_point(ts)
            points.InsertPoint(i,p)
            i=i+1
        #Insert end point of the last segment
        ts = ts_ls[-1]
        e = lgm_xml.ts_end_point(ts)
        points.InsertPoint(i,e)
        #Create spline, i.e parametric function source reresenting it
        fs = self.set_spline(points,self.resolution)
        points_array = fs.GetOutput().GetPoints()
        npoints = fs.GetOutput().GetNumberOfPoints()
        #print("NPOINTS",npoints,type(points_ls))
        #for i in range(0,npoints):
        #    print(i,points_ls.GetPoint(i))
        #Update variables for number of spline points and spline segment length
        self.set_number_of_spline_points(fs)
        self.set_spline_segment_length()
        #Create the vector of segment radiuses matching the length of spline points
        tube_radius_array = self.set_radius_data(ts_ls)
        #Color for LGAtype
        tube_lgatype_array = self.set_lgatype_data(ts_ls)
        #print("N2",fs.GetOutput().GetNumberOfPoints())
        polydata = self.set_point_polydata(fs,points_array)
        #polydata = self.set_lines_polydata(fs,lines_array)
        polydata = self.set_radius_polydata(fs,tube_radius_array,self.radius_name)
        polydata = self.set_lgatype_polydata(fs,tube_lgatype_array)
        polydata = self.set_active_scalars(polydata,self.radius_name)
        #Create tube filter (vtkTubeFilter) and add polydata into it
        #print(points_array.GetNumberOfPoints(),lines_array.GetNumberOfCells(),tube_radius_array.GetNumberOfTuples(),
        #      tube_lgatype_array.GetNumberOfTuples())
        tube_filter = self.set_tube_filter_polydata(polydata)
        #Create polydata mapper (vtkPolyDataMapper), add tube filter and set scalar ranges based on polydata
        tube_mapper = self.set_polydata_mapper_polydata(tube_filter,polydata)
        #Create actor (vtkActor) for the tube mapper
        tube_actor = self.set_tube_actor(tube_mapper)
        self.renderer.AddActor(tube_actor)
    def set_spline(self,points,resolution):
        """
        Create the spline based on the segment points in a branch
        """
        fs = vtk.vtkParametricFunctionSource()
        spline = vtk.vtkParametricSpline()
        spline.SetPoints(points)
        fs.SetParametricFunction(spline)
        fs.SetUResolution(self.resolution * points.GetNumberOfPoints())
        fs.Update()
        return fs
    def create_line(self,index1,index2):
        """
        Create lines based on the indices of the two points in the points array 
        """
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0,index1)
        line.GetPointIds().SetId(1,index2)
        return line
    def set_radius_data(self,ts_ls):
        """
        Spline points are created. Each segment is represented
        by number of spline points. The number of points depends
        on the spline resolution and segment length. For each spline
        point for the segment add segment radius
        """
        lgm_xml=LignumXML()
        tube_radius_array=vtk.vtkDoubleArray()
        tube_radius_array.SetName(self.radius_name)
        line_array = vtk.vtkCellArray()
        self.reset_index()
        for ts in ts_ls:
            l = lgm_xml.ts_length(ts)
            r = lgm_xml.ts_r(ts)
            #Calculate the number of spline points for the segment 
            nsegments=self.spline_segments(l)
            for i in range(0,nsegments):
                tube_radius_array.InsertNextValue(r)
                #index1 = self.next_index()
                #index2 = self.next_index()
                #line = self.create_line(index1,index2)
                #NOTE!!!:Why pcac83 root system needs this line_array setup
                #As local variable it will disappear. The name line_array is not used elsewhere
                #It seems when resolution (-r option) is 500 or 1000 this problem appears,
                #Update: some other roots have the same problem with differrent resolutions
                #line_array.InsertNextCell(line)
        return tube_radius_array
    def set_radius_polydata(self,fs,tube_radius_array,scalar_name):
        """
        Add the vector of radiuses to the polydata of the 
        function source (i.e. the spline)
        """
        polydata = fs.GetOutput()
        polydata.GetPointData().AddArray(tube_radius_array)
        return polydata
    def set_active_scalars(self,polydata,scalar_name):
        polydata.GetPointData().SetActiveScalars(scalar_name)
        return polydata
    def set_lgatype_data(self,ts_ls):
        lgm_xml=LignumXML()
        tube_lgatype_array=vtk.vtkDoubleArray()
        tube_lgatype_array.SetName(self.lgatype_name)
        for ts in ts_ls:
            l = lgm_xml.ts_length(ts)
            v = lgm_xml.ts_type(ts)
            #Calculate the number of spline points for the segment 
            nsegments=self.spline_segments(l)
            if v > 0:
                print('LGAvalue',v,'L',l,'N',nsegments)
            #if lgm_xml.ts_omega(ts) == 0:
            #    print(nsegments,self.spline_points,l,v,lgm_xml.ts_omega(ts))
            for i in range(0,nsegments):
                tube_lgatype_array.InsertNextValue(v)
        return tube_lgatype_array
    def set_lgatype_polydata(self,fs,tube_lgatype_array):
        polydata = fs.GetOutput()
        #polydata.GetCellData().AddArray(tube_lgatype_array)
        polydata.GetPointData().AddArray(tube_lgatype_array)
        #n = tube_lgatype_array.GetNumberOfTuples()
        #for i in range(0,n):
        #    print(i,tube_lgatype_array.GetValue(i))
        return polydata
    def set_lines_polydata(self,fs,lines_array):
        polydata = fs.GetOutput()
        polydata.SetLines(lines_array)
        return polydata
    def set_point_polydata(self,fs,points_array):
        polydata = fs.GetOutput()
        polydata.SetPoints(points_array)
        return polydata
    def set_tube_filter_polydata(self,polydata):
        tube_filter = vtk.vtkTubeFilter()
        tube_filter.AddInputData(polydata)
        tube_filter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
        tube_filter.SetNumberOfSides(self.nsides)
        tube_filter.CappingOn()
        tube_filter.Update()
        return tube_filter
    def set_polydata_mapper_polydata(self,tube_filter,tube_polydata):
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(tube_filter.GetOutputPort())
        mapper.SetScalarRange(tube_polydata.GetScalarRange())
        return mapper
    def set_tube_actor(self,tube_mapper):
        tube_actor = vtk.vtkActor()
        tube_actor.SetMapper(tube_mapper)
        return tube_actor
    def set_branch_length(self,ts_ls):
        """
        Calculate the total length of the single branch
        """
        lgm_xml=LignumXML()
        ls = [lgm_xml.ts_length(ts) for ts in ts_ls]
        length = sum(ls)
        self.branch_length=length
    def set_spline_segment_length(self):
        """
        Spline segment length is the branch length divided by number of spline points
        """
        self.spline_segment_length = self.branch_length/self.spline_points
    def spline_segments(self,segment_length):
        """
        Number of spline points for a single segments is the segment length
        divided by spline segment length rounded to the nearest integer
        """
        return int(np.round(segment_length/self.spline_segment_length))
    def set_number_of_spline_points(self,fs):
        """
        Calculate the number of spline points
        """
        n = fs.GetOutput().GetNumberOfPoints()
        self.spline_points=n
    def write_vtp_file_spline(self,file_name):
        #Reset camera is needed, otherwise the data might be not update for writer
        self.renderer.ResetCamera()
        mb_set = vtk.vtkMultiBlockDataSet()
        actors = self.renderer.GetActors()
        for (i,actor) in enumerate(actors):
            mb_set.SetBlock(i,actor.GetMapper().GetInput())
        self.writer.SetFileName(file_name)
        self.writer.SetInputData(mb_set)
        self.writer.Update()
        self.writer.Write()
        
class LignumVTKLine:
    """
    Visualize Lignum tree segments with VTK and Paraview
    Take a list of tree segments, convert to VTK representation,
    and save the output file. The idea is to create lines, 
    points and radiuses based on segment dimensions and the use vtkTubeFilter
    to create tubes around the lines (see set_polydata() and set_tube_filter()) 
    with associated scalars and colors.
    """
    def __init__(self,ts_ls,cylinder:bool,color_name:str,radius_name:str,lgatype_name:str,nsides:int):
        """
        Points, colors, lines and radiuses are side by side in their respective  arrays.
        Tricky part is to define lines and index them properly
        """
        self.ts_ls = ts_ls
        self.cylinder=cylinder
        self.nsides = nsides
        self.point_index=-1
        self.points = vtk.vtkPoints()
        #RGB colors, there is now one color array for segments but one can implement
        #as many vectors as needed. In ParaView it is then possible to visualize one
        #at a time. Please note that all vectors must have the right length with
        #values. Otherwise vtk library will crash
        self.rgb_colors = vtk.vtkUnsignedCharArray()
        self.rgb_colors.SetNumberOfComponents(3)
        self.rgb_colors.SetName(color_name)
        #Colors based on LGAtype scalar values
        self.lgatype_colors = vtk.vtkDoubleArray()
        self.lgatype_colors.SetName(lgatype_name)
        #Lines denote segment lengths and are needed for vtkTubeFiler
        self.lines = vtk.vtkCellArray()
        #For each point there is radius that vtkTubeFilter uses to calculate
        #the tree segment shape
        self.radiuses = vtk.vtkDoubleArray()
        self.radiuses.SetName(radius_name)
        #vtkPolyData will contain points and lines and their colors and radiuses
        self.lines_polydata = vtk.vtkPolyData()
        #vtkTubeFilter will display tubes based on lines and their attributes (colors, dimensions)
        self.tube_filter = vtk.vtkTubeFilter()
        #writer writes the output to ParaView compatible '.vtp' file 
        self.writer = vtk.vtkXMLPolyDataWriter()
    def get_radius_array_name(self):
        return self.radiuses.GetName()
    def get_color_array_name(self):
        return self.colors.GetName()
    def next_index(self):
        """
        Return next index for a point pointing to its place in the 
        array of points.
        """
        self.point_index = self.point_index+1
        return self.point_index
    def add_color(self,rgb_color):
        """
        Add an RGB color to the colors array
        """
        self.rgb_colors.InsertNextTypedTuple(rgb_color)
    def add_lgatype_color(self,value):
        """
        Add LGAtype scalar value
        """
        self.lgatype_colors.InsertNextTuple1(value)
    def add_point(self,point):
        """
        Add point: np.array type
        """
        self.points.InsertNextPoint(point)
    def add_radius(self,radius):
        """
        Add radius, floating point number
        """
        self.radiuses.InsertNextTuple1(radius)
    def create_line(self,index1,index2):
        """
        Create lines based on the indices of the two points in the points array 
        """
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0,index1)
        line.GetPointIds().SetId(1,index2)
        return line
    def add_line(self,line):
        self.lines.InsertNextCell(line)
    def set_polydata_points(self):
        self.lines_polydata.SetPoints(self.points)
    def set_polydata_lines(self):
        self.lines_polydata.SetLines(self.lines)
    def add_polydata_colors(self):
        #Set one scalar and the add other arrays.
        #vtkXMLPolyDataWriter then can generate vtp xml file for ParaView
        #where user can select the array to be visualized
        self.lines_polydata.GetCellData().SetScalars(self.rgb_colors)
        self.lines_polydata.GetCellData().AddArray(self.lgatype_colors)
    def add_polydata_radiuses(self):
        self.lines_polydata.GetPointData().AddArray(self.radiuses)
    def set_polydata_active_scalar(self,name):
        """
        In practice active scalar name should be tree segment radius
        """
        self.lines_polydata.GetPointData().SetActiveScalars(name)
    def set_polydata(self,active_scalar_name):
        self.set_polydata_points()
        self.set_polydata_lines()
        self.add_polydata_colors()
        self.add_polydata_radiuses()
        self.set_polydata_active_scalar(active_scalar_name)
        print(self.lines_polydata.GetNumberOfCells(),self.lgatype_colors.GetNumberOfTuples())
    def set_tube_filter(self):
        """
        Construct the tube filter, for triangulation user 
        must define number of sides, e.g. 50 seems to a good number
        """
        self.tube_filter.SetInputData(self.lines_polydata)
        self.tube_filter.SetNumberOfSides(self.nsides)
        self.tube_filter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
        self.tube_filter.CappingOn()
        self.tube_filter.Update()
    def ts_ls_to_vtp_file_line(self):
        lgm_xml = LignumXML()
        for ts in self.ts_ls:
            p = lgm_xml.ts_point(ts)
            d = lgm_xml.ts_direction(ts)
            l = lgm_xml.ts_length(ts)
            e = lgm_xml.ts_end_p(p,l,d)
            r = lgm_xml.ts_r(ts)
            rtop = lgm_xml.ts_rtop(ts)
            scalar = lgm_xml.ts_type(ts)
            if self.cylinder == True:
                rtop = r
                i1 = self.next_index()
                i2 = self.next_index()
                self.add_point(p)
                self.add_point(e)
                self.add_radius(r)
                self.add_radius(rtop)
                line = self.create_line(i1,i2)
                self.add_line(line)
                if lgm_xml.ts_omega(ts) == 0:
                    self.add_color(color_table['green'])
                else:
                    self.add_color(color_table['sandybrown'])
                self.lgatype_colors.InsertNextValue(scalar)
        self.set_polydata(self.get_radius_array_name())
        self.set_tube_filter()
    def write_vtp_file_line(self,file_name):
        """
        Write tree segments to ouput file, '.vtp' xml seems to work fine
        """
        self.ts_ls_to_vtp_file_line()
        self.writer.SetFileName(file_name)
        self.writer.AddInputDataObject(self.tube_filter.GetOutput())
        self.writer.Write()
        
class LignumXML:
    """
    Read and manipulate Lignum XML file
    """
    def __init__(self):
        pass
        #Boolean use Segment base radius as Segment top radius or not
        #self.cylinder = cylinder
        #Number of sides in segment visualization
        #self.nsides = nsides
    def ts_type(self,ts):
        """
        Return tree segment LGAtype value
        """
        children_ls = list(ts)
        child = children_ls[0]
        l = child.find('LGAtype')
        return np.float(l.text)
    def ts_omega(self,ts):
        """
        Return tree segment LGAomega value
        """
        children_ls = list(ts)
        child = children_ls[0]
        l = child.find('LGAomega')
        return np.float(l.text)
    def ts_rtop(self,ts):
        """
        Return tree segment top radius
        """
        children_ls = list(ts)
        child = children_ls[0]
        l = child.find('LGARTop')
        return np.float(l.text)
    def ts_r(self,ts):
        """
        Return tree segment radius (bottom)
        """
        children_ls = list(ts)
        child = children_ls[0]
        l = child.find('LGAR')
        return np.float(l.text)
    def ts_end_point(self,ts):
        """
        Calculate the end point
        """
        p=self.ts_point(ts)
        l=self.ts_length(ts)
        d=self.ts_direction(ts)
        return self.ts_end_p(p,l,d)
    def ts_end_p(self,p,l,d):
        """
        Calculate tree segment end point as numpy array
        """
        return np.add(p,np.multiply(l,d))
    def ts_length(self,ts):
        """
        Reutrn tree segment length
        """
        children_ls = list(ts)
        child = children_ls[0]
        l = child.find('LGAL')
        return np.float(l.text)
    def ts_direction(self,ts):
        """
        Return tree segment direction as numpy array
        """
        children_ls = list(ts)
        child = children_ls[0]
        d = child.find('direction')
        return np.fromstring(d.text,sep=' ')
    def ts_point(self,ts):
        """
        Return segment point (bottom) as numpy array
        """
        children_ls = list(ts)
        child = children_ls[0]
        p = child.find('point')
        return np.fromstring(p.text,sep=' ')
    def lignum_xml_to_treesegments(self,f):
        """
        Read Lignum XML file and return the list of all tree segments
        """
        t = ET()
        t.parse(f)
        ts_it = t.iter('TreeSegment')
        ts_ls = list(ts_it)
        return ts_ls
 
class LignumBranches:
    def __init__(self,ts_ls):
        #Dictionary of segments partitioned based on their branching order
        self.branch_dict={}
        #The list of segments in the tree
        self.ts_ls = ts_ls
    def reconstruct_branches(self):
        """
        Construct a dictionary of segments based on their branching order
        Construct a set of base points of the segments
        Find the end segments for each branch (end points are not in the point set)
        For each end segment find the preceeding segments in the same branch
        Collect the branch lists into a list and return that list of branch lists
        """
        lgm_xml=LignumXML()
        #The result will be a list of branch lists
        branch_lss=[]
        self.partition_branch_orders()
        p_set=self.construct_point_set()
        be_ls = self.find_branch_ends(p_set)
        #Construct the branch for each end segment
        for ts in be_ls:
            omega = lgm_xml.ts_omega(ts)
            b_omega_ls = self.branch_dict[omega]
            ls = self.construct_one_branch(ts,b_omega_ls)
            branch_lss.append(ls)
        return branch_lss
    def construct_one_branch(self,ts,b_omega_ls):
        """
        Recursively construct one branch (i.e. sequence of segments of a certain branching order)
        by starting from the end segment and finding backwards subsequent segments by matching 
        segment end point of the preceeding segment and the base point of its subsequent segment.
        """
        lgm_xml=LignumXML()
        ls = [ts1 for ts1 in b_omega_ls if self.point_equal(lgm_xml.ts_point(ts),lgm_xml.ts_end_point(ts1))]
        if not ls:
            #No more preceeding segments, problem solved
            return [ts]
        else:
            #Preceeding segment found
            ts2 = ls[0]
            #Remove the list element (a bit optional, but makes the list shorter anyway)
            index = b_omega_ls.index(ts2)
            b_omega_ls.pop(index)
            #Solve the problem recursively for the preceeding segment
            #and append the result to the subsequent segment  
            return self.construct_one_branch(ts2,b_omega_ls)+[ts]
    def partition_branch_orders(self):
        lgm_xml=LignumXML()
        for ts in self.ts_ls:
            self.branch_dict.setdefault(lgm_xml.ts_omega(ts),[]).append(ts)
    def find_branch_ends(self,p_set:set):
        """
        Find for each branching order the branch ends
        Return the list of terminating segments
        """
        lgm_xml=LignumXML()
        key_ls = list(self.branch_dict.keys())
        be_ls=[]
        for key in key_ls:
            #All segments of one branching order
            ls = self.branch_dict[key]
            #Terminating segments
            ls1 = [ts for ts in ls if not self.set_membership(ts,p_set)]
            #Non-terminating segments
            ls2 = [ts for ts in ls if self.set_membership(ts,p_set)]
            #Update the list of segments containing now only non-terminating segments 
            self.branch_dict[key] = ls2
            be_ls = be_ls+ls1
        return be_ls
    def construct_point_set(self):
        """
        Construct the set of base points of the segments
        Segments with end point not in the set are terminating segments
        NOTE: numpy.array is not hashable and to construct the point set
        convert numpy.array to tuple.
        """
        lgm_xml=LignumXML()
        ls = [tuple(lgm_xml.ts_point(ts)) for ts in self.ts_ls]
        return set(ls)
    def set_membership(self,ts,p_set,eps=1.0e-7):
        """
        End points and base points of segments do not match excactly
        Use eps value as the distance where the two points are considered equal
        """
        lgm_xml=LignumXML()
        dist=self.calculate_minimum_distance(lgm_xml.ts_end_point(ts),list(p_set))
        #print(dist)
        if dist <= eps:
            return True
        else:
            return False
    def point_equal(self,p1,p2,eps=1.0e-7):
        dist = distance.euclidean(p1,p2)
        if dist <= eps:
            return True
        else:
            return False
    def calculate_minimum_distance(self,p,p_ls):
        lgm_xml=LignumXML()
        ls = [(distance.euclidean(tuple(p),x)) for x in p_ls]
        return np.min(ls)
    
if __name__ == "__main__":
    parser = OP()
    parser.add_option("-i","--fxml",dest="f1",help="Read Lignum xml file")
    parser.add_option("-o","--fvtp",dest="f2",help="VTP output file")
    parser.add_option("-n","--lgatype",dest="lgatype",
                      help="Name alias for LGAtype variable, appears in color legend)")
    parser.add_option("-c","--cylinder",dest="cylinder",action="store_true",default=False,
                      help="Use segment base radius as segment top radius (pure cylinder)")
    parser.add_option("-s","--spline",dest="spline",action="store_true",default=True,
                      help="Use splines to connect segments (default)")
    parser.add_option("-l","--line",dest="line",action="store_true",default=False,
                      help="Use lines to show individual segments")
    parser.add_option("-r","--resolution",dest="resolution",help="Spline segment resolution (default=100)")
    (options,arg) = parser.parse_args()
    if options.f1 == None:
        print("No Lignum xml file")
        quit()
    if options.f2 == None:
        print("No VTP output file")
        quit()
    cylinder=False
    if options.cylinder == True:
        #Use pure cylinder, i.e. segment base radius == segment top radius
        cylinder = True
    lgatype_name = "LGAtype"
    if options.lgatype:
        lgatype_name = options.lgatype
    resolution=100
    if options.resolution:
        resolution=int(options.resolution)
    lgm_xml = LignumXML()
    ts_ls = lgm_xml.lignum_xml_to_treesegments(options.f1)
    if options.line:
        print("Using lines to show individual segments")
        lgm_line = LignumVTKLine(ts_ls,True,'SegmentColors','SegmentRadius',lgatype_name,20)
        lgm_line.write_vtp_file_line(options.f2)
    else:
        print("Using splines to connect segments")
        lgm_branches = LignumBranches(ts_ls)
        branch_lss = lgm_branches.reconstruct_branches()
        lgm_spline = LignumVTKSpline(branch_lss,True,'SegmentColors','SegmentRadius',lgatype_name,resolution,20)
        lgm_spline.tree_to_tube_actors()
        lgm_spline.write_vtp_file_spline(options.f2)
