#%% Imports
import SimpleITK as sitk
import networkx as nx
import numpy as np
from skan.csr import skeleton_to_csgraph
from scipy.stats import variation
import copy
import skimage as sk
import os

def binary_fractal_window_series(img, window_sizes):
    fractal_dim_series = []
    n_grid_samples = []
    for grid_spacing in window_sizes:
        img_arr = sitk.GetArrayFromImage(img)
        split_axis = list(img_arr.shape)
        for n,n_split in enumerate(split_axis):
            split_axis[n] =  (np.floor(n_split/grid_spacing), n_split%grid_spacing)

        img_arr = img_arr[split_axis[0][1]:, split_axis[1][1]:, split_axis[2][1]:]
        test_var = cubify(img_arr, [grid_spacing, grid_spacing, grid_spacing])

        sample_window_range = np.mean(test_var, axis=1)
        count = np.count_nonzero(sample_window_range)

        fractal_dim_series.append(count)
    return fractal_dim_series

def junction_node_subgraph(graph):
    """
    :param graph:
    :return:
    """
    jgraph = copy.deepcopy(graph)
    # Iterate through the edges and remove attributes
    for u, v, attrs in jgraph.edges(data=True):
        for attr_key in list(attrs.keys()):
            jgraph[u][v].pop(attr_key)

    parallel_edges = True
    while parallel_edges:
        parallel_edges = False
        for node in graph.nodes():
            if jgraph.has_node(node):
                if graph.degree[node] == 2: # for degree 2 nodes

                    new_edge = []
                    # assert len(list(jgraph.neighbors((item[0])))) == 2
                    for n_node in jgraph.neighbors(node):
                        new_edge.append(n_node)
                    # print(f"Processing node: {node}, with neighbors: {new_edge}")
                    if len(new_edge) == 2:
                        new_path = []
                        for edge in jgraph.edges(node):
                            if jgraph.get_edge_data(edge[0], edge[1]) != {}:
                                # print(jgraph.get_edge_data(edge[0], edge[1]))
                                new_path = jgraph.get_edge_data(edge[0], edge[1])['path']
                        jgraph.remove_node(node)
                        u, v = new_edge
                        # print(f'adding edge: {u,v}')
                        if not jgraph.has_edge(u, v):
                            jgraph.add_edge(u, v)
                            if jgraph.get_edge_data(u, v) == {}:
                                new_path.append(node)
                                nx.set_edge_attributes(jgraph, {(u, v): new_path}, name='path')
                        else:
                            parallel_edges = True
                    elif len(new_edge) == 1:
                        for edge_node in [node, new_edge[0]]:
                            if graph.degree[edge_node] == 2:
                                jgraph.remove_node(edge_node)
                #else:
                    # print(f"not processing node: {node}, with degree: {jgraph.degree[node]}, and neighbors: {list(jgraph.neighbors(node))}")

    return jgraph


def surface_area(input_im):
    """
    :param input_im: sitk.Image
    :return:
    """
    native_size = input_im.GetSize()
    mirror_image_boundary = sitk.MirrorPadImageFilter()
    mirror_image_boundary.SetPadLowerBound([1, 1, 1])
    mirror_image_boundary.SetPadUpperBound([1, 1, 1])
    padded_im = mirror_image_boundary.Execute(input_im)
    kernel = sitk.Image([3, 3, 3], sitk.sitkUInt8)
    kernel[1, 1, 0] = 1
    kernel[2, 1, 1] = 1
    kernel[0, 1, 1] = 1
    kernel[1, 2, 1] = 1
    kernel[1, 1, 2] = 1
    kernel[1, 0, 1] = 1

    conv_im = sitk.Convolution(padded_im, kernel)
    extract_filter = sitk.ExtractImageFilter()
    extract_filter.SetSize(native_size)
    extract_filter.SetIndex([1, 1, 1])

    out_im = extract_filter.Execute(conv_im)
    out_im.CopyInformation(input_im)
    # print(sitk.GetArrayFromImage(out_im).sum())
    out_im = out_im * (input_im == 0)
    area = sitk.GetArrayFromImage(out_im).flatten().sum()
    area = area * (input_im.GetSpacing()[0] * input_im.GetSpacing()[1])

    return area

def cubify(arr, newshape):
    oldshape = np.array(arr.shape)
    repeats = (oldshape / newshape).astype(int)
    tmpshape = np.column_stack([repeats, newshape]).ravel()
    order = np.arange(len(tmpshape))
    order = np.concatenate([order[::2], order[1::2]])
    # newshape must divide oldshape evenly or else ValueError will be raised
    return arr.reshape(tmpshape).transpose(order).reshape(-1, np.prod(newshape))

def Lacunarity(img, grid_spacings):
    lacunarity = []
    for grid_spacing in grid_spacings:
        img_arr = sitk.GetArrayFromImage(img)

        split_axis = list(img_arr.shape)
        for n,n_split in enumerate(split_axis):
            split_axis[n] =  (np.floor(n_split/grid_spacing), n_split%grid_spacing)

        img_arr = img_arr[split_axis[0][1]:, split_axis[1][1]:, split_axis[2][1]:]
        test_var = cubify(img_arr, [grid_spacing, grid_spacing, grid_spacing])

        test_var = np.sum(test_var, axis=1)/(grid_spacing**3)

        lacunarity.append(variation(test_var)**2)

        # print(f"The lacunarity for the length-scale of {grid_spacing}, is {variation(test_var)**2}")
    return lacunarity

def efficient_largest_ccmp_filter(image: sitk.Image):
    ccmp = sitk.ConnectedComponent(image)
    stats = sitk.LabelShapeStatisticsImageFilter()
    stats.Execute(ccmp)
    label_sizes = []
    for label in stats.GetLabels():
        label_sizes.append(stats.GetNumberOfPixels(label))
    label_max = label_sizes.index(max(label_sizes)) + 1
    image = ccmp == label_max
    return image


def get_node_degree_set(graph, degreeValue):
    deg = nx.degree(graph)
    degree_node_set = []
    for node in graph.nodes:
        if deg[node] == degreeValue:
            degree_node_set.append(node)
    return degree_node_set
# for sample, case, resolution in zip(samples, cases, resolutions):

def process_file(seg, spacing, smoothing_radius = 2, smoothing = False, downsample_radius = None, ):

    print(f"Processing with smoothing: {smoothing}, downsampling: {downsample_radius!=None}, sample spacing: {spacing}")

    grid_spacings = [2,4,8,16,32]
    results = {}
    seg = sitk.Cast(seg, sitk.sitkUInt8)
    if downsample_radius:
        downsample_radius = int(downsample_radius)
        dss_filter = sitk.ResampleImageFilter()
        dss_filter.SetInterpolator(sitk.sitkNearestNeighbor)
        dss_filter.SetSize([int(x / downsample_radius) for x in seg.GetSize()])
        dss_filter.SetOutputSpacing(
            [x * downsample_radius for x in spacing])  # sct resolution hardcoded to be 0.8 microns isotropic

        dss_filter.SetOutputOrigin(seg.GetOrigin())
        seg = dss_filter.Execute(seg)
    else:
        seg.SetSpacing(spacing)


    if smoothing:
        closing_filter = sitk.BinaryMorphologicalClosingImageFilter()
        closing_filter.SetKernelType(sitk.sitkBall)
        closing_filter.SetKernelRadius((smoothing_radius,) * seg.GetDimension())

        opening_filter = sitk.BinaryMorphologicalOpeningImageFilter()
        opening_filter.SetKernelType(sitk.sitkBall)
        opening_filter.SetKernelRadius((smoothing_radius,) * seg.GetDimension())

        seg = closing_filter.Execute(seg)
        seg = opening_filter.Execute(seg)

    seg = sitk.BinaryFillhole(seg)
    seg = efficient_largest_ccmp_filter(seg)
    spacing = seg.GetSpacing()[0]

    ccmp_villae = sitk.ConnectedComponent(seg)
    ccmp_intervillae = sitk.ConnectedComponent(seg == 0)
    results["Villous CCMP"] = sitk.GetArrayViewFromImage(ccmp_villae).max()
    results["Intervillous CCMP"] = sitk.GetArrayViewFromImage(ccmp_intervillae).max()

    results['Area'] = surface_area(seg)
    results['Volume'] = np.sum(sitk.GetArrayFromImage(seg).flatten()) * spacing ** 3
    results['Porosity'] = results['Volume'] / (seg.GetNumberOfPixels() * spacing ** 3)

    seg_arr = sitk.GetArrayFromImage(seg > 0)
    skel = sk.morphology.skeletonize_3d(seg_arr)

    # %% convert Skeleton to graph

    print('Creating graph from skeleton')
    dimensions = 3
    pixel_graph, coordinates = skeleton_to_csgraph(skel)
    coordinates = np.array(coordinates).astype(int)
    if coordinates.shape[1] > 3:
        coordinates = coordinates.T

    # coordinates = np.vstack([x for x in coordinates]).T
    # deleting the skeleton from memory, again in case we are dealing with large datasets that take up heaps of memory
    del skel

    nx_graph = nx.Graph(pixel_graph)

    # %% Terminal branch analysis
    degree_one_nodes = get_node_degree_set(nx_graph, 1)
    results['Number of Terminal Branches'] = len(degree_one_nodes)
    J_graph = junction_node_subgraph(nx_graph)
    results['Number of Branch points'] = J_graph.number_of_nodes() - results['Number of Terminal Branches']
    results['Number of Branch points'] = len(degree_one_nodes)
    results['Number of Cycles'] = len(nx.cycle_basis(J_graph))

    y = binary_fractal_window_series(ccmp_villae, grid_spacings)
    x = grid_spacings

    logx = np.log(x)
    logy = np.log(y)

    xrange = np.array(np.linspace(x[0], x[-1]))
    log_xrange = np.log(xrange)

    m, b = np.polyfit(logx, logy, 1)
    results["Fractal Dimension"] = -1*m
    # Lacunarity

    lacunarity_value = Lacunarity(ccmp_villae, grid_spacings)

    results["Lacunarity"] = sum(lacunarity_value)/len(grid_spacings)

    return results

if __name__ == '__main__':
    # FILE I/O parameters
    file_path = os.getcwd() + r'\input\67_10562_seg.nii.gz'
    spacing = [0.8,]*3

    input_img = sitk.ReadImage(file_path)
    input_img = input_img > 0
    output = process_file(input_img, spacing)
    with open(os.getcwd() + r'\output\output.txt', 'w') as f:
        for key in output.keys():
            f.write(f"{key}: {output[key]}\n")
            print(f"{key}: {output[key]}")
        print(f"Surface area to volume ratio: {output['Area']/output['Volume']}")
    print('Morphological processing complete')