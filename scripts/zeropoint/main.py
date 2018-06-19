from zeropoint import ZeroPoint
import os.path

root = "../../standards/"

if os.listdir(root) != []:
    for date in sorted(os.listdir(root)):
        if os.path.isdir(os.path.join(root, date)) and date[:3] == "201":
            zeroPoint = ZeroPoint(root, date)
            zeroPoint.ZeroPoint()
