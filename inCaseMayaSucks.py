import math
import itertools
import numpy as np
from scipy.sparse import csgraph


# referenced and modified from: https://github.com/TanaTanoi/as-rigid-as-possible-deformation
class ARAP:
    def neighboursOf(self, vertID):
        return np.where(self.neighbourMatrix[vertID] == 1)[0]

    def assignWeightForPair(self, i, j):
        if self.weightMatrix[j, i] == 0:
            weightIJ = self.weightForPair(i, j)
        else:
            weightIJ = self.weightMatrix[j, i]
        self.weightSum[i, i] += weightIJ * 0.5
        self.weightSum[j, j] += weightIJ * 0.5
        self.weightMatrix[i, j] = weightIJ

    def weightForPair(self, i, j):
        localTris = []
        for triID in self.vertsToTris[i]:
            tri = self.tris[triID]
            if i in tri and j in tri:
                localTris.append(tri)

        vertexI = self.verts[i]
        vertexJ = self.verts[j]

        cotThetaSum = 0
        for tri in localTris:
            otherVertID = list(set(tri) - set([i, j]))[0]
            otherVertex = self.verts[otherVertID]

            vA = vertexI - otherVertex
            vB = vertexJ - otherVertex
            cosTheta = vA.dot(vB) / (np.linalg.norm(vA) * np.linalg.norm(vB))
            theta = math.acos(cosTheta)

            cotThetaSum += math.cos(theta) / math.sin(theta)

        return cotThetaSum * 0.5

    def input(self, verts, tris):
        self.verts = np.array(verts)
        self.tris = np.array(tris)

        self.n = len(self.verts)
        self.vertsPrime = np.array(self.verts)
        self.vertsToTris = [[j for j, tri in enumerate(self.tris) if i in tri] for i in range(self.n)]

        self.vertsPrime = np.asmatrix(self.vertsPrime)
        self.neighbourMatrix = np.zeros((self.n, self.n))
        self.neighbourMatrix[
            tuple(zip(
                *itertools.chain(
                    *map(
                        lambda tri: itertools.permutations(tri, 2),
                        self.tris
                    )
                )
            ))
        ] = 1

        self.cellRotations = np.zeros((self.n, 3, 3))

        self.weightMatrix = np.zeros((self.n, self.n), dtype=np.float)
        self.weightSum = np.zeros((self.n, self.n), dtype=np.float)

        for vertID in range(self.n):
            neighbours = self.neighboursOf(vertID)
            for neighbourID in neighbours:
                self.assignWeightForPair(vertID, neighbourID)

    def update(self, fixedIDs, handleIDs, deformationMatrices):
        deformationMatrices = list(map(np.matrix, deformationMatrices))
        self.deformationVerts = []
        for i in range(self.n):
            if i in handleIDs:
                deformedVector = np.append(self.verts[i], 1)
                deformedVector = deformedVector.dot(deformationMatrices[handleIDs.index(i)])
                deformedVector = np.delete(deformedVector, 3).flatten()
                deformedVector = np.squeeze(np.asarray(deformedVector))
                self.deformationVerts.append((i, deformedVector))
            elif i in fixedIDs:
                self.deformationVerts.append((i, self.verts[i]))

        # extended laplacian matrix
        deformationVertsNum = len(self.deformationVerts)
        self.laplacianMatrix = np.zeros([self.n + deformationVertsNum] * 2, dtype=np.float32)
        self.laplacianMatrix[:self.n, :self.n] = csgraph.laplacian(self.weightMatrix)
        for i in range(deformationVertsNum):
            vertID = self.deformationVerts[i][0]
            ni = i + self.n
            self.laplacianMatrix[ni, vertID] = 1
            self.laplacianMatrix[vertID, ni] = 1

        # precompute PiArray
        self.PiArray = []
        for i in range(self.n):
            vertI = self.verts[i]
            neighbourIDs = self.neighboursOf(i)
            neighboursNum = len(neighbourIDs)

            Pi = np.zeros((3, neighboursNum))

            for ni in range(neighboursNum):
                nID = neighbourIDs[ni]

                vertJ = self.verts[nID]
                Pi[:, ni] = (vertI - vertJ)
            self.PiArray.append(Pi)

    def calculateCellRotations(self):
        for vertID in range(self.n):
            rotation = self.calculateRotationMatrixForCell(vertID)
            self.cellRotations[vertID] = rotation

    def applyCellRotations(self):
        for i in range(self.n):
            self.bArray[i] = np.zeros((1, 3))
            neighbours = self.neighboursOf(i)
            for j in neighbours:
                wij = self.weightMatrix[i, j] / 2.0
                rij = self.cellRotations[i] + self.cellRotations[j]
                pij = self.verts[i] - self.verts[j]
                self.bArray[i] += (wij * rij.dot(pij))

        self.vertsPrime = np.linalg.solve(self.laplacianMatrix, self.bArray)[:self.n]

    def calculateRotationMatrixForCell(self, vertID):
        covarianceMatrix = self.calculateConvarianceMatrixForCell(vertID)

        U, s, VTranspose = np.linalg.svd(covarianceMatrix)

        rotation = VTranspose.T.dot(U.T)
        if np.linalg.det(rotation) <= 0:
            U[:0] *= -1
            rotation = VTranspose.T.dot(U.T)
        return rotation

    def calculateConvarianceMatrixForCell(self, vertID):
        vertIPrime = self.vertsPrime[vertID]

        neighbourIDs = self.neighboursOf(vertID)
        neighboursNum = len(neighbourIDs)

        Di = np.zeros((neighboursNum, neighboursNum))

        Pi = self.PiArray[vertID]
        PiPrime = np.zeros((3, neighboursNum))

        for ni in range(neighboursNum):
            nID = neighbourIDs[ni]

            Di[ni, ni] = self.weightMatrix[vertID, nID]

            vertJPrime = self.vertsPrime[nID]
            PiPrime[:, ni] = (vertIPrime - vertJPrime)

        PiPrime = PiPrime.T

        return Pi.dot(Di).dot(PiPrime)

    def apply(self, iterations):
        deformationVertsNum = len(self.deformationVerts)

        self.bArray = np.zeros((self.n + deformationVertsNum, 3))
        for i in range(deformationVertsNum):
            self.bArray[self.n + i] = self.deformationVerts[i][1]

        for t in range(iterations):
            self.calculateCellRotations()
            self.applyCellRotations()


def testAllTheTests(verts, tris, fixedIDs, handleIDs, inputMatrices, iterations):
    arap = ARAP()
    arap.input(
        verts=verts,  # vertex positions as np array with shape (-1, 3)
        # vertex indices per triangle as np array with shape (-1, 3)
        tris=np.array(tris).reshape(-1, 3)
    )
    arap.update(
        fixedIDs=fixedIDs,  # vertex IDs that shouldn't move
        handleIDs=handleIDs,  # vertex IDs that you want to move
        # targets transforms per handle (vertexes you want to move in order)
        deformationMatrices=inputMatrices,
    )
    arap.apply(iterations)
    return arap.vertsPrime  # returns the resulting verex positions


if __name__ == '__main__':
    print(testAllTheTests(
        # 4 vertex positions (x, y, z)
        [
            [-1., .0, .0],
            [1., .0, .0],
            [.0, .0, -1.],
            [.0, .0, 1.]
        ],
        [[0, 1, 2], [1, 2, 3]],  # triangles from vertex IDs (zero indexed)
        [],  # no fixed verts
        [0],  # first vert is being moved
        # identity matrix isn't very meaningful here I guess, but a good test
        np.identity(4).reshape(1, 4, 4),
        1  # 1 iteration
    ))  # this should print the same vertex positions with some precision loss
