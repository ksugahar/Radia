import os

import win32com.client


source_pptx = os.environ.get("RADIA_MCP_PPTX_EXPORT_SOURCE")
if not source_pptx:
    raise SystemExit("Set RADIA_MCP_PPTX_EXPORT_SOURCE to a local .pptx file.")

PPT = win32com.client.Dispatch("PowerPoint.Application")
PPT.Visible = 1
PPT.Presentations.Open(source_pptx)
PPT.ActivePresentation.Slides(1).Export(os.getcwd() + "/slide_001.png", "png")
PPT.ActivePresentation.Slides(1).Export(os.getcwd() + "/slide_001.jpg", "jpg")
PPT.ActivePresentation.Close()
PPT.Quit()
