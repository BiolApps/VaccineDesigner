import re
import pandas as pd


def blast_parsing_fnc(res_path):
    with open(res_path, "r") as file:
        blastp_output = file.read()
    # Regular expression patterns
    header_pattern = re.compile(r"Sequences producing significant alignments:\s+\(Bits\)\s+Value\n([\s\S]+?)\n\n", re.DOTALL)
    #alignment_pattern = re.compile(r"^(\S+)\s+RecName: Full=([^[]+)\[.*?\]\nLength=(\d+).*?\nScore =\s+(\S+).*?Expect =\s+(\S+)", re.MULTILINE)
    # Find header and alignments
    header_match = header_pattern.search(blastp_output)
    header_text = header_match.group(1)
    info=header_text.split('\n')
    info=[i for i in info if i!='']
    queries=[i.split(' ')[0] for i in info]
    queries_text=['>'+ i for i in queries]


    data={'Query':[],'Rec-name':[],'Flags':[],'Short':[],'Length':[],'Score':[],'Expect':[],'identities':[],'Positives':[],'Gaps':[]}
    for query in queries_text:
        index1=blastp_output.find(query)
        index2=blastp_output[index1:].find('Query')
        query_string=blastp_output[index1:index1+index2]



        # Define regular expression patterns to extract information
        rec_name_pattern = re.compile(r"RecName: Full=(.*?)(?=;|\n)")
        flags_pattern = re.compile(r"Flags: (.*?)(?=Length=|\Z)", re.DOTALL)
        length_pattern = re.compile(r"Length=(\d+)")
        score_pattern = re.compile(r"Score = (.+?) bits")
        expect_pattern = re.compile(r"Expect = (.+?),")
        identities_pattern = re.compile(r"Identities = (.+?),")
        positives_pattern = re.compile(r"Positives = (.+?),")
        gaps_pattern = re.compile(r"Gaps = (.+?)\n")
        shorts_pattern=re.compile(r"Short=(.*?);")

        rec_name_match = rec_name_pattern.search(query_string)
        flags_match = flags_pattern.search(query_string)
        length_match = length_pattern.search(query_string)
        score_match = score_pattern.search(query_string)
        expect_match = expect_pattern.search(query_string)
        identities_match = identities_pattern.search(query_string)
        positives_match = positives_pattern.search(query_string)
        gaps_match = gaps_pattern.search(query_string)
        short_match=shorts_pattern.search(query_string)
        data['Query'].append(query[1:])
        
        if rec_name_match:
            rec_name = rec_name_match.group(1)
            data['Rec-name'].append(rec_name)
        else:
            data['Rec-name'].append('')

        if flags_match:
            flags = flags_match.group(1)
            data['Flags'].append(flags)
        else:
            data['Flags'].append('')

        if length_match:
            length = length_match.group(1)
            data['Length'].append(length)
        else:
            data['Length'].append('')

        if score_match:
            score = score_match.group(1)
            data['Score'].append(score)
        else:
            data['Score'].append('')

        if expect_match:
            expect = expect_match.group(1)
            data['Expect'].append(expect)
        else:
            data['Expect'].append('')

        if identities_match:
            identities = identities_match.group(1)
            data['identities'].append(identities)
        else:
            data['identities'].append('')

        if positives_match:
            positives = positives_match.group(1)
            data['Positives'].append(positives)
        else:
            data['Positives'].append('')

        if gaps_match:
            gaps = gaps_match.group(1)
            data['Gaps'].append(gaps)
        else:
            data['Gaps'].append('')

        if short_match:
            short=short_match.group(1)
            data['Short'].append(short)
        else:
            data['Short'].append('')
    df=pd.DataFrame.from_dict(data)
    return df
  
  

def blast_parsing_fnc_meta(res_path):
    with open(res_path, "r") as file:
        blastp_output = file.read()
    # Regular expression patterns
    header_pattern = re.compile(r"Sequences producing significant alignments:\s+\(Bits\)\s+Value\n([\s\S]+?)\n\n", re.DOTALL)
    #alignment_pattern = re.compile(r"^(\S+)\s+RecName: Full=([^[]+)\[.*?\]\nLength=(\d+).*?\nScore =\s+(\S+).*?Expect =\s+(\S+)", re.MULTILINE)
    # Find header and alignments
    header_match = header_pattern.search(blastp_output)
    header_text = header_match.group(1)
    info=header_text.split('\n')
    info=[i for i in info if i!='']
    queries=[i.split(' ')[0] for i in info]
    queries_text=['>'+ i for i in queries]


    data={'Query':[],'Rec-name':[],'Length':[],'Score':[],'Expect':[],'identities':[],'Positives':[],'Gaps':[]}
    for query in queries_text:
        index1=blastp_output.find(query)
        index2=blastp_output[index1:].find('Query')
        query_string=blastp_output[index1:index1+index2]



        # Define regular expression patterns to extract information
        rec_name_pattern = re.compile(r"{qr} (.*?)(?=;|\n)".format(qr=query))
        
        length_pattern = re.compile(r"Length=(\d+)")
        score_pattern = re.compile(r"Score = (.+?) bits")
        expect_pattern = re.compile(r"Expect = (.+?),")
        identities_pattern = re.compile(r"Identities = (.+?),")
        positives_pattern = re.compile(r"Positives = (.+?),")
        gaps_pattern = re.compile(r"Gaps = (.+?)\n")
        

        rec_name_match = rec_name_pattern.search(query_string)
        
        length_match = length_pattern.search(query_string)
        score_match = score_pattern.search(query_string)
        expect_match = expect_pattern.search(query_string)
        identities_match = identities_pattern.search(query_string)
        positives_match = positives_pattern.search(query_string)
        gaps_match = gaps_pattern.search(query_string)
        
        data['Query'].append(query[1:])
        
        if rec_name_match:
            rec_name = rec_name_match.group(1)
            data['Rec-name'].append(rec_name)
        else:
            data['Rec-name'].append('')

        

        if length_match:
            length = length_match.group(1)
            data['Length'].append(length)
        else:
            data['Length'].append('')

        if score_match:
            score = score_match.group(1)
            data['Score'].append(score)
        else:
            data['Score'].append('')

        if expect_match:
            expect = expect_match.group(1)
            data['Expect'].append(expect)
        else:
            data['Expect'].append('')

        if identities_match:
            identities = identities_match.group(1)
            data['identities'].append(identities)
        else:
            data['identities'].append('')

        if positives_match:
            positives = positives_match.group(1)
            data['Positives'].append(positives)
        else:
            data['Positives'].append('')

        if gaps_match:
            gaps = gaps_match.group(1)
            data['Gaps'].append(gaps)
        else:
            data['Gaps'].append('')

        
    df=pd.DataFrame.from_dict(data)
    return df
    
