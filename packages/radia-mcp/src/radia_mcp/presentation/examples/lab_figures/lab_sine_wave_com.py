import glob
import os
import win32com.client
from numpy import *

PPT = win32com.client.Dispatch("PowerPoint.Application")
PPT.Visible = 1
PPT.Presentations.Add();
PPT.ActivePresentation.PageSetup.SlideSize = 7
PPT.ActivePresentation.PageSetup.SlideWidth = 2000
PPT.ActivePresentation.PageSetup.SlideHeight = 500

PPT.ActivePresentation.PageSetup.FirstSlideNumber = 1
layout = PPT.ActivePresentation.SlideMaster.CustomLayouts.Item(7)
newSlide = PPT.ActivePresentation.Slides.AddSlide(1,layout)

myFreeForm = PPT.ActiveWindow.Selection.SlideRange.Shapes.BuildFreeform(0,50,200)
x = arange(50,450,1)
y = -100*sin( 4*pi*(x-x[0])/200.)+200

for n in range(len(x)):
	myFreeForm.AddNodes(0,0,x[n],y[n])
myFreeForm.ConvertToShape()

PPT.ActivePresentation.SaveAs(os.getcwd() + u"test.pptx")

