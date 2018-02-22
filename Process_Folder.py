# @File(label = "Input directory", style = "directory") srcFile
# @File(label = "Output directory", style = "directory") dstFile
# @String(label = "File extension", value=".tif") ext
# @String(label = "File name contains", value = "") containString
# @boolean(label = "Keep directory structure when saving", value = true) keepDirectories

# See also Process_Folder.ijm for a version of this code
# in the ImageJ 1.x macro language.

import os
from ij import IJ, ImagePlus
from ij.plugin.frame import RoiManager
from ij.gui import Roi



def run():
  srcDir = srcFile.getAbsolutePath()
  dstDir = dstFile.getAbsolutePath()
  for root, directories, filenames in os.walk(srcDir):
    filenames.sort();
    for filename in filenames:
      # Check for file extension
      if not filename.endswith(ext):
        continue
      # Check for file name pattern
      if containString not in filename:
        continue
      process(srcDir, dstDir, root, filename, keepDirectories)
 
def process(srcDir, dstDir, currentDir, fileName, keepDirectories):
  print "Processing:"
  ROIpath,V1S1RL = os.path.split(currentDir)
  # Opening the image
  print "Open image file", fileName
  imp = IJ.openImage(os.path.join(currentDir, fileName))
  imp.show()
  print "Compute F/F0"
  IJ.run(imp,"F div F0","6")
  imp.close()
  
  V1S1RL = "RoiSet " + V1S1RL + ".zip"
  # Put your processing commands here!
  ROIPath = ROIpath + "/" +V1S1RL
  
  print "Open ROIS", ROIPath
  rm = RoiManager.getInstance()
  if not rm:
    rm = RoiManager()
  rm.reset()
  rm.runCommand("open",ROIPath)
  rm.runCommand("Show All")
  print "Multi Measure"
  imp = IJ.getImage()
  rt = rm.multiMeasure(imp)
  imp.close()
  imp = IJ.getImage()
  imp.close()

  # Saving the image
  saveDir = currentDir.replace(srcDir, dstDir) if keepDirectories else dstDir
  if not os.path.exists(saveDir):
    os.makedirs(saveDir)
  print "Saving to", saveDir
  rt.save(os.path.join(saveDir, fileName + '.csv'));
  # IJ.saveAs(imp, "Tiff", os.path.join(saveDir, fileName));
  imp.close()
 
run()
