import sys
import math
import itertools
import json
import numpy as np
from scipy.sparse import csgraph
import maya.OpenMayaMPx as OpenMayaMPx
import maya.OpenMaya as OpenMaya
import maya.cmds as MC


pluginNodeName = 'arapDeformer'
pluginNodeID = OpenMaya.MTypeId(0x77113355)

# reference: http://discourse.techart.online/t/maya-2016-python-deformer-plugin/5239/6
kApiVersion = MC.about(apiVersion=True)
if kApiVersion < 201600:
    kInput = OpenMayaMPx.cvar.MPxDeformerNode_input
    kInputGeom = OpenMayaMPx.cvar.MPxDeformerNode_inputGeom
    kOutputGeom = OpenMayaMPx.cvar.MPxDeformerNode_outputGeom
    kEnvelope = OpenMayaMPx.cvar.MPxDeformerNode_envelope
else:
    kInput = OpenMayaMPx.cvar.MPxGeometryFilter_input
    kInputGeom = OpenMayaMPx.cvar.MPxGeometryFilter_inputGeom
    kOutputGeom = OpenMayaMPx.cvar.MPxGeometryFilter_outputGeom
    kEnvelope = OpenMayaMPx.cvar.MPxGeometryFilter_envelope


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
            zip(
                *itertools.chain(
                    *map(
                        lambda tri: itertools.permutations(tri, 2),
                        self.tris
                    )
                )
            )
        ] = 1

        self.cellRotations = np.zeros((self.n, 3, 3))

        self.weightMatrix = np.zeros((self.n, self.n), dtype=np.float)
        self.weightSum = np.zeros((self.n, self.n), dtype=np.float)

        for vertID in range(self.n):
            neighbours = self.neighboursOf(vertID)
            for neighbourID in neighbours:
                self.assignWeightForPair(vertID, neighbourID)

    def update(self, fixedIDs, handleIDs, deformationMatrices):
        deformationMatrices = map(np.matrix, deformationMatrices)
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


class ArapDeformerNode(OpenMayaMPx.MPxDeformerNode):
    # Amount to push the vertices by
    attrFixedIDsStr = OpenMaya.MObject()
    attrHandleIDsStr = OpenMaya.MObject()
    attrInputMatrices = OpenMaya.MObject()
    attrIterations = OpenMaya.MObject()

    def __init__(self):
        OpenMayaMPx.MPxDeformerNode.__init__(self)
        self.once = True

    def postConstructor(self):
        selfObj = self.thisMObject()
        selfFn = OpenMaya.MFnDependencyNode(selfObj)
        self.inputMatricesPlug = selfFn.findPlug('inputMatrices')
        self.fixedIDsStrPlug = selfFn.findPlug('fixedIds')
        self.handleIDsStrPlug = selfFn.findPlug('handleIds')

    def deform(self, data, geomIter, localToWorldMatrix, geomIDx):
        # envelope = data.inputValue(kEnvelope).asFloat()
        if self.once:
            inputGeomObj = self.getInputGeom(data, geomIDx)
            counts = OpenMaya.MIntArray()
            tris = OpenMaya.MIntArray()
            verts = OpenMaya.MPointArray()
            mesh = OpenMaya.MFnMesh(inputGeomObj)
            mesh.getTriangles(counts, tris)
            mesh.getPoints(verts, OpenMaya.MSpace.kObject)

            self.arapd = ARAP()
            self.arapd.input(
                verts=[[verts[i].x, verts[i].y, verts[i].z] for i in range(verts.length())],
                tris=np.array(tris).reshape(-1, 3)
            )
            self.once = False

        dataStr = self.fixedIDsStrPlug.asString()
        try:
            fixedIDs = json.loads(dataStr)
        except Exception:
            fixedIDs = None

        dataStr = self.handleIDsStrPlug.asString()
        try:
            handleIDs = json.loads(dataStr)
        except Exception:
            handleIDs = None

        inputMatrices = []
        for i in xrange(self.inputMatricesPlug.numElements()):
            matPlug = self.inputMatricesPlug.elementByPhysicalIndex(i)
            if matPlug.isConnected():
                matrixObj = matPlug.asMObject()
                matrixData = OpenMaya.MFnMatrixData(matrixObj)
                matrix = matrixData.matrix()
                inputMatrices.append(
                    [
                        [matrix(0, j) for j in range(4)],
                        [matrix(1, j) for j in range(4)],
                        [matrix(2, j) for j in range(4)],
                        [matrix(3, j) for j in range(4)],
                    ]
                )

        if inputMatrices and handleIDs and len(handleIDs) == len(inputMatrices):
            if not fixedIDs:
                fixedIDs = []
            self.arapd.update(
                fixedIDs=fixedIDs,
                handleIDs=handleIDs,
                deformationMatrices=inputMatrices,
            )
            iterationHandle = data.inputValue(ArapDeformerNode.attrIterations)
            iterations = iterationHandle.asInt()
            self.arapd.apply(iterations)

            while not geomIter.isDone():
                idx = geomIter.index()
                geomIter.setPosition(OpenMaya.MPoint(*self.arapd.vertsPrime[idx]))
                geomIter.next()

    def getInputGeom(self, data, geomIDx):
        inputHandle = data.outputArrayValue(kInput)
        inputHandle.jumpToElement(geomIDx)
        inputGeomObj = inputHandle.outputValue().child(kInputGeom).asMesh()
        return inputGeomObj


def nodeCreator():
    return OpenMayaMPx.asMPxPtr(ArapDeformerNode())


def nodeInitializer():
    numAttrFn = OpenMaya.MFnNumericAttribute()

    # Setup attributes
    ArapDeformerNode.attrIterations = numAttrFn.create(
        'iterations',
        'iterations',
        OpenMaya.MFnNumericData.kInt,
        1
    )
    numAttrFn.setMin(1)
    numAttrFn.setMax(100)
    numAttrFn.setChannelBox(True)
    ArapDeformerNode.addAttribute(ArapDeformerNode.attrIterations)

    matrixAttrFn = OpenMaya.MFnMatrixAttribute()
    typedAttrFn = OpenMaya.MFnTypedAttribute()

    strFn = OpenMaya.MFnStringData()
    defaultValue = strFn.create('null')
    ArapDeformerNode.attrFixedIDsStr = typedAttrFn.create('fixedIds', 'fixedIds', OpenMaya.MFnData.kString, defaultValue)
    typedAttrFn.setReadable(True)
    typedAttrFn.setWritable(True)
    typedAttrFn.setStorable(True)
    ArapDeformerNode.addAttribute(ArapDeformerNode.attrFixedIDsStr)

    strFn = OpenMaya.MFnStringData()
    defaultValue = strFn.create('null')
    ArapDeformerNode.attrHandleIDsStr = typedAttrFn.create('handleIds', 'handleIds', OpenMaya.MFnData.kString, defaultValue)
    typedAttrFn.setReadable(True)
    typedAttrFn.setWritable(True)
    typedAttrFn.setStorable(True)
    ArapDeformerNode.addAttribute(ArapDeformerNode.attrHandleIDsStr)

    ArapDeformerNode.attrInputMatrices = matrixAttrFn.create('inputMatrices', 'inputMatrices')
    matrixAttrFn.setArray(True)
    matrixAttrFn.setReadable(True)
    matrixAttrFn.setWritable(True)
    matrixAttrFn.setStorable(True)
    ArapDeformerNode.addAttribute(ArapDeformerNode.attrInputMatrices)

    # Link inputs that change the output of the mesh
    ArapDeformerNode.attributeAffects(
        ArapDeformerNode.attrIterations,
        kOutputGeom
    )
    ArapDeformerNode.attributeAffects(
        ArapDeformerNode.attrFixedIDsStr,
        kOutputGeom
    )
    ArapDeformerNode.attributeAffects(
        ArapDeformerNode.attrHandleIDsStr,
        kOutputGeom
    )
    ArapDeformerNode.attributeAffects(
        ArapDeformerNode.attrInputMatrices,
        kOutputGeom
    )


def initializePlugin(mobject):
    mplugin = OpenMayaMPx.MFnPlugin(mobject)
    try:
        mplugin.registerNode(
            pluginNodeName,
            pluginNodeID,
            nodeCreator,
            nodeInitializer,
            OpenMayaMPx.MPxNode.kDeformerNode
        )
    except Exception:
        sys.stderr.write('Failed to register node: ' + pluginNodeName)
        raise


def uninitializePlugin(mobject):
    mplugin = OpenMayaMPx.MFnPlugin(mobject)
    try:
        mplugin.deregisterNode(pluginNodeID)
    except Exception:
        sys.stderr.write('Failed to deregister node: ' + pluginNodeName)
        raise
