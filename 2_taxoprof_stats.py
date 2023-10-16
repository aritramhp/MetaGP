import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse as ap

import util


def merging_abundance(indir):
    tax_file = os.path.join(indir,'OTUtable.rel_abundance.tab')
    cmd_merge = 'merge_metaphlan_tables.py ' + indir + '/profiles/*.txt > ' + tax_file
    os.system(cmd_merge)
    return tax_file


def separate_taxrank(merged_file, indir):
    outdir = os.path.join(indir,'Taxonomic_binning')
    util.create_dir(outdir)
    
    fp_k = open(outdir+'/0_kingdom.tab', 'w')
    fp_p = open(outdir+'/1_phylum.tab', 'w')
    fp_c = open(outdir+'/2_class.tab', 'w')
    fp_o = open(outdir+'/3_order.tab', 'w')
    fp_f = open(outdir+'/4_family.tab', 'w')
    fp_g = open(outdir+'/5_genera.tab', 'w')
    fp_s = open(outdir+'/6_species.tab', 'w')
    fp_t = open(outdir+'/7_strain.tab', 'w')

    fp = open(merged_file,'r')
    line = fp.readline()
    line = fp.readline()
    fp_k.write(line)
    fp_p.write(line)
    fp_c.write(line)
    fp_o.write(line)
    fp_f.write(line)
    fp_g.write(line)
    fp_s.write(line)
    fp_t.write(line)

    while True:
        line = fp.readline()
        if line:
            # tax = line.split('\t')[0].split('|')[-1]
            tax = line.split('|')[-1]
            if 'k__' in tax:
                fp_k.write(tax)
            elif 'p__' in tax:
                fp_p.write(tax)
            elif 'c__' in tax:
                fp_c.write(tax)
            elif 'o__' in tax:
                fp_o.write(tax)
            elif 'f__' in tax:
                fp_f.write(tax)
            elif 'g__' in tax:
                fp_g.write(tax)
            elif 's__' in tax:
                fp_s.write(tax)
            elif 't__' in tax:
                fp_t.write(tax)
        else:
            fp.close()
            fp_k.close()
            fp_p.close()
            fp_c.close()
            fp_o.close()
            fp_f.close()
            fp_g.close()
            fp_s.close()
            fp_t.close()
            break
    return outdir

def plot_relabundance(abun_file):
    colors269 = ["#000000","#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#997D87","#5A0007","#809693","#FEFFE6","#1B4400","#4FC601","#3B5DFF","#4A3B53","#FF2F80","#61615A","#BA0900","#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA","#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299","#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2","#C2FF99","#001E09","#00489C","#6F0062","#0CBD66","#EEC3FF","#456D75","#B77B68","#7A87A1","#788D66","#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459","#456648","#0086ED","#886F4C","#34362D","#B4A8BD","#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F","#938A81","#575329","#00FECF","#B05B6F","#8CD0FF","#3B9700","#04F757","#C8A1A1","#1E6E00","#7900D7","#A77500","#6367A9","#A05837","#6B002C","#772600","#D790FF","#9B9700","#549E79","#FFF69F","#201625","#72418F","#BC23FF","#99ADC0","#3A2465","#922329","#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98","#A4E804","#324E72","#6A3A4C","#83AB58","#001C1E","#D1F7CE","#004B28","#C8D0F6","#A3A489","#806C66","#222800","#BF5650","#E83000","#66796D","#DA007C","#FF1A59","#8ADBB4","#1E0200","#5B4E51","#C895C5","#320033","#FF6832","#66E1D3","#CFCDAC","#D0AC94","#7ED379","#012C58","#7A7BFF","#D68E01","#353339","#78AFA1","#FEB2C6","#75797C","#837393","#943A4D","#B5F4FF","#D2DCD5","#9556BD","#6A714A","#001325","#02525F","#0AA3F7","#E98176","#DBD5DD","#5EBCD1","#3D4F44","#7E6405","#02684E","#962B75","#8D8546","#9695C5","#E773CE","#D86A78","#3E89BE","#CA834E","#518A87","#5B113C","#55813B","#E704C4","#00005F","#A97399","#4B8160","#59738A","#FF5DA7","#F7C9BF","#643127","#513A01","#6B94AA","#51A058","#A45B02","#1D1702","#E20027","#E7AB63","#4C6001","#9C6966","#64547B","#97979E","#006A66","#391406","#F4D749","#0045D2","#006C31","#DDB6D0","#7C6571","#9FB2A4","#00D891","#15A08A","#BC65E9","#FFFFFE","#C6DC99","#203B3C","#671190","#6B3A64","#F5E1FF","#FFA0F2","#CCAA35","#374527","#8BB400","#797868","#C6005A","#3B000A","#C86240","#29607C","#402334","#7D5A44","#CCB87C","#B88183","#AA5199","#B5D6C3","#A38469","#9F94F0","#A74571","#B894A6","#71BB8C","#00B433","#789EC9","#6D80BA","#953F00","#5EFF03","#E4FFFC","#1BE177","#BCB1E5","#76912F","#003109","#0060CD","#D20096","#895563","#29201D","#5B3213","#A76F42","#89412E","#1A3A2A","#494B5A","#A88C85","#F4ABAA","#A3F3AB","#00C6C8","#EA8B66","#958A9F","#BDC9D2","#9FA064","#BE4700","#658188","#83A485","#453C23","#47675D","#3A3F00","#061203","#DFFB71","#868E7E","#98D058","#6C8F7D","#D7BFC2","#3C3E6E","#D83D66","#2F5D9B","#6C5E46","#D25B88","#5B656C","#00B57F","#545C46","#866097","#365D25","#252F99","#00CCFF","#674E60","#FC009C","#92896B"]
    colors17= ['#3eb489','#ff6ec7','#ffd12b','#03324a','#16621c','#5f2e4c','#f9584b','#596fff','#81dd4d','#e6b710','#ff8980','#99adb7','#b5db52','#3fe0d2','#c90076','#2ffd51','#ff8933']
    colors = ['#556b2f','#8b4513','#483d8b','#008000','#00008b','#daa520','#008b8b','#7f007f','#8fbc8f','#b03060','#ffff00',
              '#00ff00','#8a2be2','#ff0000','#00ff7f','#00ffff','#0000ff','#adff2f','#dc143c','#da70d6','#ff7f50','#ff00ff',
              '#f0e68c','#6495ed','#90ee90','#ff1493','#87cefa','#fff8dc','#ffb6c1']
    abun_tab = pd.read_table(abun_file, sep='\t')
    abun_tab = abun_tab.set_index('clade_name')
    new_df = abun_tab[abun_tab.apply(lambda row: all(float(column) <= 1 for column in row), axis=1)]
    df1 = abun_tab.drop(new_df.index, errors='ignore')
    df_others = new_df.sum(axis='index')
    df_others.name = 'Others'
    df1 = df1.append(df_others)
    df1 = df1.T
    no_of_color = df1.shape[1]
    if no_of_color-1 <= len(colors):
        clrs = colors[:no_of_color-1]
    else:
        clrs = colors + colors[:no_of_color-1-len(colors)]
    df1.plot(kind="bar",stacked = True, color = clrs+['#aeaeae'],figsize=(100,30),position=.5)
    plt.xticks(fontsize=30, rotation=45)
    plt.yticks(fontsize=30)
    plt.legend(bbox_to_anchor=(0.5, -0.2), fontsize="30", loc='upper center', borderaxespad=0, ncol=10)
    plt.savefig(abun_file.replace('.tab','.png'))

   
parser = ap.ArgumentParser()
parser.add_argument('-i','--indir',dest='input_dir', type=str, required=True, help='Input directory')

args = parser.parse_args()
indir = args.input_dir
indir = os.path.join(indir,'2_taxonomic_profile')
tax_file = merging_abundance(indir)
tax_bining = separate_taxrank(tax_file, indir)
for f in os.listdir(tax_bining):
    if f.endswith('.tab'):
        filename = os.path.join(tax_bining,f)
        plot_relabundance(filename)