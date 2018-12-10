import json
import os.path


def GetSettings(userfile='Analysis/config/user.json'):
    with open('Analysis/config/defaults.json') as F:
        settings = json.load(F)

    if not os.path.isfile(userfile):
        return settings

    with open(userfile) as F:
        userSettings = json.load(F)

    for uSection in userSettings:
        if uSection not in settings:
            settings[uSection] = userSettings[uSection]
            continue

        for uKey in userSettings[uSection]:
            settings[uSection][uKey] = userSettings[uSection][uKey]

    return settings
