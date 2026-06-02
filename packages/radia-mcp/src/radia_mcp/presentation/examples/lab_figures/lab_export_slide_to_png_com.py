import glob
import os
import win32com.client

PPT = win32com.client.Dispatch("PowerPoint.Application")
PPT.Visible = 1
PPT.Presentations.Open(u"w:/2019_テーマと委員.pptx")
PPT.ActivePresentation.Slides(1).Export(os.getcwd() + u"/ksugahar.png","png");
PPT.ActivePresentation.Slides(1).Export(os.getcwd() + u"/ksugahar.jpg","jpg");
PPT.ActivePresentation.Close
PPT.Quit

