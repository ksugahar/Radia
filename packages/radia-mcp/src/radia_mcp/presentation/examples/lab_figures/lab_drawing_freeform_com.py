import glob
import os
import win32com.client
from numpy import *

PPT = win32com.client.Dispatch("PowerPoint.Application")
PPT.Visible = 1
PPT.Presentations.Add();
PPT.ActivePresentation.PageSetup.SlideSize = 7;
PPT.ActivePresentation.PageSetup.SlideWidth = 2000;
PPT.ActivePresentation.PageSetup.SlideHeight = 500;

PPT.ActivePresentation.PageSetup.FirstSlideNumber = 1;
layout = PPT.ActivePresentation.SlideMaster.CustomLayouts.Item(7);
newSlide = PPT.ActivePresentation.Slides.AddSlide(1,layout);

for n in range(1,99):
	a = PPT.ActiveWindow.Selection.SlideRange.Shapes.AddShape(n,30*n-10,300,30,30)
	a = PPT.ActiveWindow.Selection.SlideRange.Shapes.AddShape(14,30*n-10,330,60,30)
	a.Fill.ForeColor.RGB = 255*255*255+3*255*255+3*255
	a.TextFrame.TextRange.Text = str(n)
	a.TextFrame.TextRange.TextColor = 0

myFreeForm = PPT.ActiveWindow.Selection.SlideRange.Shapes.BuildFreeform(0,50,200);

x = linspace(50,450,401)
y = 50*sin(4*pi*(x-x[0])/100)+200

for n in range(len(x)):
	myFreeForm.AddNodes(0,0,x[n],y[n])

myFreeForm.ConvertToShape()

PPT.ActivePresentation.SaveAs(os.getcwd() + "utest.pptx")

