import sys
import os
import urllib.parse
import urllib.request

url = 'http://www.uniprot.org/uniprot/'

values = {
'format':'gff',
'query':'lalala'
}

data = urllib.parse.urlencode(values)
data = data.encode('utf8')
req = urllib.request.Request(url, data)
response = urllib.request.urlopen(req)
the_page = response.read()

print(the_page)



#url_values = urllib.parse.urlencode(values)
#full_url = url + '?' + url_values
#data = urllib.request.open(full_url)

#data = urlparse(params)
#request = urllib.request(url, data)
#contact = "" # Please set your email address here to help us debug in case of problems.
#request.add_header('User-Agent', 'Python %s' % contact)
#response = urllib.urlopen(request)
#page = response.read()
