import json
import maya.mel as MM
import maya.cmds as MC
from PySide import QtGui


def __load_plugin():
    pluginPath = 'arapDeformerNode.py'
    if not MC.pluginInfo(pluginPath, query=True, loaded=True):
        if MC.loadPlugin(pluginPath) is None:
            MC.warning('Failed to load plugin.')
            return False
    return True


def __apply(mesh):
    name = mesh.split('|')[-1] + '_arapDeformer'
    name = MM.eval('formValidObjectName("%s");' % name)
    return MC.deformer(mesh, name=name, type='arapDeformer')


def applyArapDeformer(mesh):
    if not __load_plugin():
        return

    return __apply(mesh)


class ArapDeformerDialog(QtGui.QDialog):

    instance = None

    def __init__(self, *args):
        '''
        build UI
        '''
        super(ArapDeformerDialog, self).__init__(*args)

        self.setWindowTitle('ARAP Deformer Tool')
        self.setMinimumSize(128, 128)

        mainLayout = QtGui.QVBoxLayout()
        self.setLayout(mainLayout)

        addHandleButton = QtGui.QPushButton('Add Handle')
        mainLayout.addWidget(addHandleButton)
        addHandleButton.clicked.connect(self.addHandleCB)

    def addHandleCB(self, *args):
        selection = MC.ls(selection=True)
        selectedVertIDs = []
        meshPath = None
        for vtxSel in selection:
            if meshPath is None:
                meshPath = vtxSel.split('.')[0]
            if meshPath != vtxSel.split('.')[0]:
                MC.warning('select vertices from one mesh!')
                return
            if '[' in vtxSel:
                if ':' in vtxSel:
                    start, end = vtxSel.split('[')[-1].split(']')[0]. split(':')
                    start = int(start)
                    end = int(end)
                    selectedVertIDs.extend(range(start, end + 1))
                else:
                    selectedVertIDs.append(int(vtxSel.split('[')[-1].split(']')[0]))
            else:
                MC.warning('Select vertices only!')
                return
            if vtxSel.split('.')[-1].split('[')[0] != 'vtx':
                MC.warning('Select vertices only!')
                return

        bbox = MC.exactWorldBoundingBox(selection)
        centerPos = (
            (bbox[0] + bbox[3]) / 2.,
            (bbox[1] + bbox[4]) / 2.,
            (bbox[2] + bbox[5]) / 2.
        )
        handleControl = MC.spaceLocator(position=(0, 0, 0), name='handleTemp')[0]
        MC.move(centerPos[0], centerPos[1], centerPos[2], handleControl)
        MC.makeIdentity(handleControl, apply=True, t=True, r=True, s=True, normal=False)

        arapDeformer = MC.listConnections('%s.inMesh' % meshPath, type='arapDeformer')
        if arapDeformer:
            arapDeformer = arapDeformer[0]
        else:
            arapDeformer = applyArapDeformer(meshPath)[0]

        try:
            handleIDs = json.loads(MC.getAttr('%s.handleIds' % arapDeformer))
        except Exception:
            handleIDs = []
        if handleIDs is None:
            handleIDs = []
        existingHandleCount = len(handleIDs)

        for i in range(len(selectedVertIDs)):
            MC.connectAttr('%s.worldMatrix[0]' % handleControl, '%s.inputMatrices[%s]' % (arapDeformer, i + existingHandleCount), force=True)

        handleIDs.extend(selectedVertIDs)
        MC.setAttr('%s.handleIds' % arapDeformer, json.dumps(handleIDs), type='string')


def showTool():
    if ArapDeformerDialog.instance is None:
        ArapDeformerDialog.instance = ArapDeformerDialog()
    ArapDeformerDialog.instance.show()
