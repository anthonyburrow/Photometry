from __future__ import print_function
from googleapiclient.discovery import build
from httplib2 import Http
from oauth2client import file, client, tools

# Get access to API
SCOPES = 'https://www.googleapis.com/auth/drive.metadata.readonly'
store = file.Storage('etc/token.json')
creds = store.get()
if not creds or creds.invalid:
    flow = client.flow_from_clientsecrets('etc/client_id.json', SCOPES)
    creds = tools.run_flow(flow, store)
service = build('drive', 'v3', http=creds.authorize(Http()))


def getfolderid(path):
    folder_query = "name = '%s' and mimeType = '%s'" % (path,
                                                        "application/vnd.google-apps.folder")
    results = service.files().list(q=folder_query).execute()
    files = results.get('files', [])
    for f in files:
        if f['name'] == path:
            return f['id']


print(getfolderid('Photometry'))
