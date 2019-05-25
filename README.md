as-rigid-as-possible (arap) maya deformer, in python

![keps.000.gif](https://github.com/leventt/arap/blob/master/keps/keps.000.gif "keps.000.gif")


Intention of this implementation was first to understand how it works by prototyping it in python and then share it as a reference, if someone may find it useful, as I do, to read it in python.
It is important to mention that inCaseMayaSucks.py is pretty much a rewrite of this:
https://github.com/TanaTanoi/as-rigid-as-possible-deformation
So I get a sense for while re-writing it.

Having said that, there are way more performant implementations out there and I would hate if you lost time trying to make this one work to find out it doesn't perform as fast for dense meshes as you might expect.

I would actually recommend the libigl implementation here:
https://github.com/libigl/libigl/blob/master/include/igl/arap.h
or Shizuo Kaji has some great implementations, if you want to use something similar and compare:
https://github.com/shizuo-kaji/CageDeformerMaya
https://github.com/shizuo-kaji/ProbeDeformerMaya
https://github.com/shizuo-kaji/PoissonMLS
(last one isn't ARAP but MLS)

You may appreciate this demo as well:
https://github.com/libigl/libigl/blob/master/tutorial/406_FastAutomaticSkinningTransformations/main.cpp

Overall, I felt that MLS (Moving Least Squares) method was more robust and faster.

If you don't have maya, check out inCaseMayaSucks.py
That should work with python3, scipy and numpy, in which case miniconda is great to run that:
https://docs.conda.io/en/latest/miniconda.html

inCaseMayaSucks.py also has some adjustments (converting itertool objects to lists and tuples) to make it work in python3.



That out of the way...
If you want to try this one out and you have a maya installation...
1.) You would need to install numpy and scipy for maya's python interpreter. This is challenging if you are on Windows, so here is a link, before you may waste time trying to figure that out:
https://forums.autodesk.com/t5/maya-programming/guide-how-to-install-numpy-scipy-in-maya-windows-64-bit/td-p/5796722
(There are much cleaner ways to install this I feel like but this may be straight-forward)
2.) You would need to put the plug-ins folder into your maya plugins path. Straight forward way to do this is to find the maya folder on your system either in your home folder or Documents on Windows... Here are more ways to do it: https://knowledge.autodesk.com/support/maya/learn-explore/caas/CloudHelp/cloudhelp/2016/ENU/Maya/files/GUID-FA51BD26-86F3-4F41-9486-2C3CF52B9E17-htm.html
3.) You will need to also put arapDeformer.py into a maya scripts path... Most straight-forward one is that same maya folder mentioned before/scripts...
(MAYA_SCRIPT_PATH mentioned here: https://knowledge.autodesk.com/support/maya/learn-explore/caas/CloudHelp/cloudhelp/2015/ENU/Maya/files/Environment-Variables-File-path-variables-htm.html)
4.) After you can successfully import that arapDeformer module... call the showTool function.



Also check this prototype if you may want to use the libigl library from python, like I was trying before for a prototype:
https://github.com/leventt/elastik

Cheers \o/
