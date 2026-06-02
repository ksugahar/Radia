import glob
import os
import win32com.client
from numpy import *

os.system('taskkill.exe /f /im POWERPNT.EXE')

objPPT = win32com.client.Dispatch("PowerPoint.Application")
objPPT.Visible = 1
objPPT.Presentations.Add();
objPPT.ActivePresentation.PageSetup.SlideSize = 7;
objPPT.ActivePresentation.PageSetup.SlideWidth  = 80/25.4*72;
objPPT.ActivePresentation.PageSetup.SlideHeight = 80/25.4*72;

objPPT.ActivePresentation.PageSetup.FirstSlideNumber = 1;
layout = objPPT.ActivePresentation.SlideMaster.CustomLayouts.Item(7);
newSlide = objPPT.ActivePresentation.Slides.AddSlide(1,layout);

N = 9
M = 9

for n in range(0,N):
	for m in range(0,M):
		x0 = (80*n/N)/25.4*72;
		y0 = (80*m/M)/25.4*72;
		lx = (80/N)/25.4*72;
		ly = (80/M)/25.4*72;
		objRect = objPPT.ActiveWindow.Selection.SlideRange.Shapes.AddShape(1,x0,y0,lx,ly)
		objRect.Fill.ForeColor.RGB = int(255*n/N)*256*256 + 0*256*256 + int(255*m/M)
		objRect.Line.Weight = 0
		objText = objRect.TextFrame.TextRange.InsertBefore(str(n+1) + 'x' + str(m+1))
		objText.Characters().Font.Size = 6

objPPT.ActivePresentation.Slides(1).Export(os.getcwd() + u"/table.png","png",800,800);



#		objRect.Fill.ForeColor.RGB =   0*256*256 +   0*256 + 255  #red
#		objRect.Fill.ForeColor.RGB =   0*256*256 + 255*256 +   0  #green
#		objRect.Fill.ForeColor.RGB = 255*256*256 +   0*256 +   0  #blue

#for n in range(1,99):
#	a = objPPT.ActiveWindow.Selection.SlideRange.Shapes.AddShape(n,30*n-10,300,30,30)
#	a = objPPT.ActiveWindow.Selection.SlideRange.Shapes.AddShape(14,30*n-10,330,60,30)
#	a.Fill.ForeColor.RGB = 255*255*255+3*255*255+3*255
#	a.TextFrame.TextRange.Text = str(n)
#	a.TextFrame.TextRange.TextColor = 0
#
#myFreeForm = objPPT.ActiveWindow.Selection.SlideRange.Shapes.BuildFreeform(0,50,200);
#
#x = linspace(50,450,401)
#y = 50*sin(4*pi*(x-x[0])/100)+200
#  for n in range(len(x)):
#	myFreeForm.AddNodes(0,0,x[n],y[n])
#
#myFreeForm.ConvertToShape()
#
#objPPT.ActivePresentation.SaveAs(os.getcwd() + r"test.pptx")

