

class CenterLineDetermination:
    def __init__(self, mandrel_name=None):

        import trimesh
        import matplotlib.pyplot as plt
        from trimesh.path.simplify import resample_spline
        import numpy as np

        if mandrel_name is None:
            import easygui
            path = easygui.fileopenbox()

        self.mandrel_name = mandrel_name

        path = mandrel_name + '.stl'

        mesh = trimesh.load(path)

        self.mesh = mesh

        mesh = mesh.process()

        sliced_mesh = mesh.section(plane_origin=[0, 0, 0], plane_normal=[0, 1, 0])

        sliced_mesh2_d = sliced_mesh.to_planar()[0]

        HTM = sliced_mesh.to_planar()[1]

        line = []
        nb_line = len(sliced_mesh2_d.entities)
        for i in range(nb_line):
            line.append(sliced_mesh2_d.entities[i].points)
        # [5] Affichage des lignes
        plt.figure()
        for i in range(nb_line):
            plt.plot(sliced_mesh2_d.vertices[line[i][:], 0], sliced_mesh2_d.vertices[line[i][:], 1], label=i,
                     marker='+')
        plt.legend()
        plt.show(block=True)

        groupe1 = input('Numeros du 1er groupe de lignes separés par une virgule sans espace,'
                        ' +/- selon le sens de parcours :')
        groupe2 = input('Numeros du 2nd groupe de lignes separés par une virgule sans espace,'
                        ' +/- selon le sens de parcours :')
        groupe1 = groupe1.split(',')
        groupe1 = [int(i) for i in groupe1]
        groupe2 = groupe2.split(',')
        groupe2 = [int(i) for i in groupe2]

        # [6] Reconstruction des 2 groupes de lignes
        if groupe1[0] >= 0:
            line1 = line[groupe1[0]][:]
        else:
            line1 = line[-groupe1[0]][::-1]
        if len(groupe1) > 1:
            for i in range(1, len(groupe1)):
                if groupe1[i] >= 0:
                    line1 = np.hstack((line1, line[groupe1[i]][:]))
                else:
                    line1 = np.hstack((line1, line[-groupe1[i]][::-1]))

        if groupe2[0] >= 0:
            line2 = line[groupe2[0]][:]
        else:
            line2 = line[-groupe2[0]][::-1]
        if len(groupe2) > 1:
            for i in range(1, len(groupe2)):
                if groupe2[i] >= 0:
                    line2 = np.hstack((line2, line[groupe2[i]][:]))
                else:
                    line2 = np.hstack((line2, line[-groupe2[i]][::-1]))

        # [7] Suppression des elements doublés et interpolation
        line1 = [ii for n, ii in enumerate(line1) if ii not in line1[:n]]
        line2 = [ii for n, ii in enumerate(line2) if ii not in line2[:n]]
        line1 = np.array((sliced_mesh2_d.vertices[line1[:], 0], sliced_mesh2_d.vertices[line1[:], 1])).transpose()
        line1 = resample_spline(line1, count=1000, degree=1)
        line2 = np.array((sliced_mesh2_d.vertices[line2[:], 0], sliced_mesh2_d.vertices[line2[:], 1])).transpose()
        line2 = resample_spline(line2, count=1000, degree=1)

        # [8] Création de la ligne des centres
        coord = []
        for i in range(line1.shape[0]):
            coord.append((line1[i] + line2[i]) / 2)
        coord = np.asarray(coord)

        # [9] Affichage dans un graphique
        plt.figure()
        plt.plot(line1[:, 0], line1[:, 1], label='1')
        plt.plot(line2[:, 0], line2[:, - 1], label='2')
        plt.legend()
        plt.plot(coord[:, 0], coord[:, 1], label='Ligne des centres')
        # [10] Passage en 3D
        liste_centroid = np.dot(HTM, np.hstack((coord, np.zeros(coord.shape[0]).reshape((coord.shape[0], 1)),
                                                np.ones(coord.shape[0]).reshape(
                                                    (coord.shape[0], 1)))).transpose()).transpose()
        # [11] Inversion du sens de parcours de la ligne des centres
        if np.linalg.norm(liste_centroid[0, 0:3]) > np.linalg.norm(liste_centroid[-1, 0:3]):
            liste_centroid = liste_centroid[::-1, 0:3]

        self.coordinates = liste_centroid

    def export_2_txt(self):
        import numpy as np
        f = open(self.mandrel_name + '_Points.txt', 'w')

        for i in range(np.shape(self.coordinates)[0]):
            f.write('Point.' + str(i) + ';' + str(round(self.coordinates[i, 0], 6)).replace('.', ',') +
                    ';0,000000;' + str(round(self.coordinates[i, 2], 6)).replace('.', ',') + '\n')
        f.close()

    def plot3d(self):
        from mayavi import mlab

        mlab.triangular_mesh(self.mesh.vertices[:, 0],
                             self.mesh.vertices[:, 1],
                             self.mesh.vertices[:, 2],
                             self.mesh.faces,
                             color=(1, 1, 1)
                             )
        mlab.plot3d(self.coordinates[:, 0], self.coordinates[:, 1], self.coordinates[:, 2])
        mlab.show()



