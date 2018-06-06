import os
import re
import Bio
from Bio import SeqIO
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='auto blast name giver')
    parser.add_argument("-a", "--align", type=str, help='path to input file with alignments', default=os.curdir+'/short_res.txt')
    parser.add_argument("-o", "--output_file", type=str, help='path and name of output file',default=os.curdir + '/len_parsed_out.txt')
    parser.add_argument("-i", "--input_file", type=str, help='path to input file with contigs', default=os.curdir+'/pacbio-only-smrt.fasta')
    args = parser.parse_args()
    inp = args.input_file
    out = args.output_file
    aln = args.align
    parse_id = re.compile(r'''
                        Query=\s
                        (?P<id>[\d\w|]*)
                        \nfull
                        ''', re.VERBOSE)
    parse_other = re.compile(r'''
                        \-[0-9]+\-[\d\w.]*[\s\t]+(PREDICTED:\s)*
                        (Sequence\sof\s[\w]*\s[\w\d]*\sfrom\s)*
                        (Genomic\s[sS]equence\s[Ff]or\s)*
                        (?P<name>[\w]+\s[\w]+)
                        \s
                        [\s\S]*?
                        [\s\t]+
                        (?P<e_value>[\d.\-e]+)
                        [\t\s]+\n
                              ''', re.VERBOSE)
    ''''''
    with open(aln, 'r') as inp :
        text = inp.readlines()

    lst = {}

    for i in range(0,len(text)-11,12):
        new_list_seqs = []
        t = []
        temp_string = text[i] + text[i+6] + text[i+7] + text[i+8] + text[i+9] + text[i+10]
        parsed_id = parse_id.finditer(temp_string)
        parsed_other = parse_other.finditer(temp_string)
        for j in parsed_id:
            #print(j.group('id'))
            id = j.group('id')
        for y in parsed_other:
            if y.group('name') == 'Ovis canadensis':
                print(y.group('name'),y.group('e_value'))
            new_list_seqs.append([y.group('e_value'), y.group('name')])
        if len(new_list_seqs) == 0:
            continue
        first_name = new_list_seqs[0][1]
        new_list_seqs = sorted(new_list_seqs)
        temp_list = []
        temp_list.append(new_list_seqs[0])
        for i in range(len(new_list_seqs) - 1):
            if (float(new_list_seqs[i + 1][0]) - float(new_list_seqs[i][0])) > 0.0001:
                break
            else:
                temp_list.append(new_list_seqs[i + 1])
        for x in temp_list:
            t.append(x[1])
        names_set = set(t)
        names_set = list(names_set)
        counted = []
        for z in range(len(names_set)):
            x = t.count(names_set[z])
            #print('XXX',x)
            names_set[z] = [x, names_set[z]]
            counted.append(x)
        maximum = max(names_set)
        names_set.remove(maximum)
        names_set.append([0,''])
        for j in range(len(names_set)):
            if maximum[0] == names_set[j][0]:
                name = first_name
                break
            else:
                #print('MAX', maximum)
                name = maximum[1]
        #print('Species of '+ id +' is:', name)
        lst.update({id:name})
    print(lst)

    enemies_set={}

    not_interesting = ['Arabidopsis', 'Arabis', 'Boechera', 'Brassica', 'Camelina',
                                'Capsella','Heuchera','Bunias','Byrsonima',
                                'Barbarea','Nicotiana', 'Zea','Citrullus',
                                'Draba','Spinacia','Erysimum','Hibiscus',
                                'Tarenaya', 'Ricinus','Chenopodium','Lepidium',
                                'Rosa','Turritis', 'Cynara', 'Momordica',
                                'Ananas', 'Medicago', 'Sclerotinia', 'Herrania',
                                'Arachis','Coffea', 'Capsicum', 'Oryza',
                                'Helianthus','Armoracia', 'Pinus','Ipomoea',
                                'Lupinus', 'Cucumis', 'Jatropha', 'Juglans',
                                'Olea','Silene', 'Asparagus', 'Picea',
                                'Beta', 'Sisymbrium','Grias', 'Theobroma', 'Ziziphus',
                                'Pugionium', 'Descurainia', 'Petunia', 'Vigna',
                                'Sesamum', 'Citrus', 'Manihot', 'Linum', 'Triticum',
                                'Castilleja', 'Lobularia', 'Gossypium', 'Megacarpaea',
                                'Cicer', 'Populus', 'Byrsonima', 'Sinapis', 'Saccharum',
                                'Geranium', 'Daucus', 'Prunus', 'Ajuga', 'Cynomorium',
                                'Sorghum','Setaria',
                                'Cardamine','Crucihimalaya','Eutrema',
                                'Glycine', 'Lotus', 'Olimarabidopsis', 'Pachycladon',
                                'Raphanus', 'Schrenkiella', 'Solanum', 'Spirodela',
                                'Thellungiella', 'Vitis']
    not_interesting = set(not_interesting)
    interesting = {'unknown':0}

    values = []
    for value in lst.values():
        values.append(value)
    values = set(values)
    values = list(values)
    for i in values:

        genus = i.split(' ')
        #print(i)
        if genus[0] in not_interesting:
            not_interesting = list(not_interesting)
            not_interesting.append(i)
            not_interesting = set(not_interesting)
        else:
            continue
    print('HELLO', not_interesting)
    with open(inp,'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            try:
                genus = lst[record.id]
                if genus in list(not_interesting):
                    continue
                else:
                    if genus in interesting.keys():
                        interesting[genus] += len(record.seq)
                    else:
                        interesting.update({genus:len(record.seq)})
            except KeyError:
                interesting['unknown'] += len(record.seq)
                #print(record.id,' ' ,len(record.seq))

    with open(out, 'w') as out:
        total = 0
        for key,value in interesting.items():
            out.write(str(key) +' : ' + str(value)+'\n')
            total+=value
        out.write('Total : '+str(total)+'\n')
        total = total - interesting['unknown']
        out.write('Total without unknowns: ' + str(total) + '\n')