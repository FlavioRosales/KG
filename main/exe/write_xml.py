def write_xml(path, namefile_h5, dt = 0.0):
    import h5py
    import numpy as np
    import os

    os.chdir(path)

    hf = h5py.File(namefile_h5, 'r')

    name_groups = list(hf.keys())

    groups = []

    for name_group in name_groups:
        groups.append(hf.get(name_group))

    name_datasets = []

    for group in groups:
        name_datasets.append(list(group.keys()))

    Xcoord = 'Xcoord'
    Ycoord = 'Ycoord'
    Zcoord = 'Zcoord'

    print(path, namefile_h5)

    datasets = []
    for name in name_datasets:
        name.remove(Xcoord)
        name.remove(Ycoord)
        name.remove(Zcoord)
        datasets.append(name)

    name_datasets = datasets[0]

    n = len(name_datasets)
    #print(n)
    base_name = name_datasets[0].split('_')[0]

    #
    # star xmf file
    #
    namefile_xmf = namefile_h5.split(".")[0] + ".xmf"

    new_file = open(namefile_xmf,"w")

    new_file.write('<?xml version="1.0" ?>\n')
    new_file.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
    new_file.write('<Xdmf Version="2.0">\n')
    new_file.write('    <Domain>\n')
    new_file.write('        <Grid Name="CellTime" GridType="Collection" CollectionType="Temporal">\n')

    t = 0.0
    for i in range(n):
        name_dataset = base_name + '_' + str(i)
        new_file.write('            <Grid Name="Concentric refinement" GridType="Collection">\n')
        new_file.write('                    <Time Type="Single" Value="'+str(t)+'"/>\n')
        for name_group in name_groups:

            X = np.array(hf.get(name_group + '/' + Xcoord))
            Y = np.array(hf.get(name_group + '/' + Ycoord))
            Z = np.array(hf.get(name_group + '/' + Zcoord))

            dim = X.shape
            dimensions = str(dim[0]) + " " + str(dim[1]) + " " + str(dim[2])

            X = X[0,0,:]
            Y = Y[0,:,0]
            Z = Z[:,0,0]

            x = str(X).replace('[', ' ').replace(']',' ')
            y = str(Y).replace('[', ' ').replace(']',' ')
            z = str(Z).replace('[', ' ').replace(']',' ')

            new_file.write('                <Grid Name="'+name_group+'" GridType="Uniform">\n')
            new_file.write('                    <Topology TopologyType="3DRectMesh" Dimensions="'+ dimensions +'"/>\n')
            new_file.write('                    <Geometry GeometryType="VXVYVZ">\n')
            #
            # X coordinate
            #
            new_file.write('                        <DataItem Name="Xcoord" Dimensions="'+ str(X.shape[0]) +'" NumberType="Float" Precision="8" Format="XML">\n')
            new_file.write(x)
            new_file.write('                        </DataItem>\n')
            #
            # Y coordinate
            #
            new_file.write('                        <DataItem Name="Ycoord" Dimensions="'+ str(Y.shape[0]) +'" NumberType="Float" Precision="8" Format="XML">\n')
            new_file.write(y)
            new_file.write('                        </DataItem>\n')
            #
            # Z coordinate
            #
            new_file.write('                        <DataItem Name="Zcoord" Dimensions="'+ str(Z.shape[0]) +'" NumberType="Float" Precision="8" Format="XML">\n')
            new_file.write(z)
            new_file.write('                        </DataItem>\n')
            new_file.write('                    </Geometry>\n')
            #
            # Scalar
            #
            new_file.write('                <Attribute Name="' + namefile_h5.split(".")[0] + '" AttributeType="Scalar" Center="Node">\n')
            new_file.write('                    <DataItem Dimensions="'+ dimensions +'" NumberType="Float" Precision="8" Format="HDF">\n')
            new_file.write('                        '+namefile_h5 + ':/' + name_group + '/' + name_dataset + '\n')
            new_file.write('                    </DataItem>\n')
            new_file.write('                </Attribute>\n')
            new_file.write('            </Grid>\n')

            #print(namefile_h5 + '/' + name_group + '/' + name_dataset)
        #
        # next time step
        #
        new_file.write('            </Grid>')
        t = t + dt
    new_file.write('        </Grid>\n')
    new_file.write('    </Domain>\n')
    new_file.write('</Xdmf>\n')
    new_file.close()

def write_vector_xml(path, namefile_h5, dt = 0.0):
    import h5py
    import numpy as np
    import os

    os.chdir(path)

    hf = h5py.File(namefile_h5, 'r')

    name_groups = list(hf.keys())

    groups = []

    for name_group in name_groups:
        groups.append(hf.get(name_group))

    name_datasets = []

    for group in groups:
        name_datasets.append(list(group.keys()))

    Xcoord = 'Xcoord'
    Ycoord = 'Ycoord'
    Zcoord = 'Zcoord'

    datasets = []
    for name in name_datasets:
        name.remove(Xcoord)
        name.remove(Ycoord)
        name.remove(Zcoord)
        datasets.append(name)

    name_datasets = datasets[0]

    n = int(len(name_datasets) / 3)

    base_name_x = name_datasets[0].split('_')[0]
    base_name_y = name_datasets[n].split('_')[0]
    base_name_z = name_datasets[2*n].split('_')[0]

    #
    # star xmf file
    #
    namefile_xmf = namefile_h5.split(".")[0] + ".xmf"

    new_file = open(namefile_xmf,"w")

    new_file.write('<?xml version="1.0" ?>\n')
    new_file.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
    new_file.write('<Xdmf Version="2.0">\n')
    new_file.write('    <Domain>\n')
    new_file.write('        <Grid Name="CellTime" GridType="Collection" CollectionType="Temporal">\n')

    t = 0.0
    for i in range(n):

        name_dataset_x = base_name_x + '_' + str(i)
        name_dataset_y = base_name_y + '_' + str(i)
        name_dataset_z = base_name_z + '_' + str(i)

        new_file.write('            <Grid Name="Concentric refinement" GridType="Collection">\n')
        new_file.write('                    <Time Type="Single" Value="'+str(t)+'"/>\n')

        for name_group in name_groups:

            X = np.array(hf.get(name_group + '/' + Xcoord))
            Y = np.array(hf.get(name_group + '/' + Ycoord))
            Z = np.array(hf.get(name_group + '/' + Zcoord))

            dim = X.shape
            dimensions = str(dim[0]) + " " + str(dim[1]) + " " + str(dim[2])

            X = X[0,0,:]
            Y = Y[0,:,0]
            Z = Z[:,0,0]

            x = str(X).replace('[', ' ').replace(']',' ')
            y = str(Y).replace('[', ' ').replace(']',' ')
            z = str(Z).replace('[', ' ').replace(']',' ')

            new_file.write('                <Grid Name="'+name_group+'" GridType="Uniform">\n')
            new_file.write('                    <Topology TopologyType="3DRectMesh" Dimensions="'+ dimensions +'"/>\n')
            new_file.write('                    <Geometry GeometryType="VXVYVZ">\n')
            #
            # X coordinate
            #
            new_file.write('                        <DataItem Name="Xcoord" Dimensions="'+ str(X.shape[0]) +'" NumberType="Float" Precision="8" Format="XML">\n')
            new_file.write(x)
            new_file.write('                        </DataItem>\n')
            #
            # Y coordinate
            #
            new_file.write('                        <DataItem Name="Ycoord" Dimensions="'+ str(Y.shape[0]) +'" NumberType="Float" Precision="8" Format="XML">\n')
            new_file.write(y)
            new_file.write('                        </DataItem>\n')
            #
            # Z coordinate
            #
            new_file.write('                        <DataItem Name="Zcoord" Dimensions="'+ str(Z.shape[0]) +'" NumberType="Float" Precision="8" Format="XML">\n')
            new_file.write(z)
            new_file.write('                        </DataItem>\n')
            new_file.write('                    </Geometry>\n')
            #
            # vector field
            #
            print(namefile_h5, name_group, name_dataset_x)
            print(namefile_h5, name_group, name_dataset_y)
            print(namefile_h5, name_group, name_dataset_z)
            new_file.write('                <Attribute Name="' + namefile_h5.split(".")[0] + '" AttributeType="Vector" Center="Node">\n')
            new_file.write('                    <DataItem Dimensions="'+ dimensions +' 3" Function = "JOIN($0,$1,$2)" ItemType = "Function">\n')
            new_file.write('                    <DataItem Dimensions="'+ dimensions +'" NumberType="Float" Precision="8" Format="HDF">\n')
            new_file.write('                        '+namefile_h5 + ':/' + name_group + '/' + name_dataset_x + '\n')
            new_file.write('                    </DataItem>\n')
            new_file.write('                    <DataItem Dimensions="'+ dimensions +'" NumberType="Float" Precision="8" Format="HDF">\n')
            new_file.write('                        '+namefile_h5 + ':/' + name_group + '/' + name_dataset_y + '\n')
            new_file.write('                    </DataItem>\n')
            new_file.write('                    <DataItem Dimensions="'+ dimensions +'" NumberType="Float" Precision="8" Format="HDF">\n')
            new_file.write('                        '+namefile_h5 + ':/' + name_group + '/' + name_dataset_z + '\n')
            new_file.write('                    </DataItem>\n')
            new_file.write('                    </DataItem>\n')
            new_file.write('                </Attribute>\n')
            new_file.write('            </Grid>\n')
        #
        # next time step
        #
        new_file.write('            </Grid>')
        t = t + dt
    new_file.write('        </Grid>\n')
    new_file.write('    </Domain>\n')
    new_file.write('</Xdmf>\n')
    new_file.close()


path = '/home/flavio/Codes/KG/main/exe/test/'
dt = 5*6.60E-03

print('Hola mundo')

write_xml(path = path, namefile_h5 = 'phi.h5', dt = dt)
write_xml(path = path, namefile_h5 = 'phi2.h5', dt = dt)

print('Adios mundo')
