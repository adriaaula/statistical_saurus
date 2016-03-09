import re
import urllib.parse
import urllib.request
from urllib.error import  URLError

url = 'http://www.uniprot.org/uniprot/'

values = {
'format':'gff',
'query':'A0A009ESS6'
}

data = urllib.parse.urlencode(values)
data = data.encode('utf8')
req = urllib.request.Request(url, data)
try:
    response = urllib.request.urlopen(req)
except URLError as e:
    if hasattr(e, 'reason'):
        print('We failed to reach a server.')
        print('Reason: ', e.reason)
    elif hasattr(e, 'code'):
        print('The server couldn\'t fulfill the request.')
        print('Error code: ', e.code)

the_page = response.read()

the_page = str(the_page)

regex = re.compile("Transmembrane")




#url_values = urllib.parse.urlencode(values)
#full_url = url + '?' + url_values
#data = urllib.request.open(full_url)

#data = urlparse(params)
#request = urllib.request(url, data)
#contact = "" # Please set your email address here to help us debug in case of problems.
#request.add_header('User-Agent', 'Python %s' % contact)
#response = urllib.urlopen(request)
#page = response.read()
