


import requests

membrane = "MADSNGTITVEELKKLLEQWNLVIGFLFLTWICLLQFAYANRNRFLYIIKLIFLWLLWPVTLACFVLAAVYRINWITGGIAIAMACLVGLMWLSYFIASFRLFARTRSMWSFNPETNILLNVPLHGTILTRPLLESELVIGAVILRGHLRIAGHHLGRCDIKDLPKEITVATSRTLSYYKLGASQRVAGDSGFAAYSRYRIGNYKLNTDHSSSSDNIALLVQ"




def checkMHCi(sequence,method='recommended',allele='HLA-A*01:01',length=11):
    data_pack = {'method':method, 'sequence_text':sequence, 'allele':allele, 'length':length}

    #POST request the API and collect the response in Response obj
    response = requests.post('http://tools-cluster-interface.iedb.org/tools_api/mhci/',data=data_pack)

    #print(response.text)

    # Split along new lines
    content = response.text.split('\n')

    return content

print(checkMHCi("HTTDPSFLGRY"))
