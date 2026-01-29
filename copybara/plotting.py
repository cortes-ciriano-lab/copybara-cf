"""
Script to plot final copy number output
Created: 20/02/2025
Python 3.9.7
Carolin Sauer
"""

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 5})
import random

import copybara.cn_functions as cnfitter

# plotting
def plotting(bin_data, seg_data, offsets, cat_colours, meta, outdir):
# def plotting(bin_data, seg_data, offsets, goi_list, cat_colours, gene_colour, meta, outdir):
    # Prepare figure
    title = f"{meta['sample']} | fitted purity={meta['fitted_purity']}; fitted ploidy={meta['ploidy']}; PAG={meta['pag']}"
    subtitle = f"CN deviation={meta['purity']}; error function={meta['gof']}; coverage={meta['coverage']}x"
    fig_name = f"{outdir}/{meta['sample']}_copy_number_plot.pdf"
    plt.figure(figsize=(3.5, 1.25))
    # prepare chromosome input fpr plotting
    chromosomes = list(dict.fromkeys(offsets))
    chr_label = [x.replace('chr','') for x in chromosomes]
    midpoints = []
    # chr_limits = [offsets['chr1']['start']]
    chr_limits = []
    # # gois
    # gois_plotting = []
    # for r in goi_list:
    #     gchr,gene = r[0],r[3]
    #     gpos = int(r[2])-int(r[1])+1 + offsets[gchr]['offset']
    #     gois_plotting.append([gpos, gene])
    # Plot each chromosome
    for chrom in chromosomes:
        # offsets
        x_offset = offsets[chrom]['offset']
        midpoints.append(offsets[chrom]['mid'])
        chr_limits.append(offsets[chrom]['end'])            
        # Filter data for the current chromosome
        bins = [d for d in bin_data if d['chrom'] == chrom]
        segments = [d for d in seg_data if d['chrom'] == chrom]
        # Plot Binned Data
        for cat, col in cat_colours.items():
            subset = [d for d in bins if d['category'] == cat]
            POS = [x_offset + d['pos'] for d in subset]
            CN = [d['binned'] for d in subset]
            # scatter bins
            plt.scatter(
                POS,
                CN,
                color=col,
                s=1,
                alpha=0.5
            )
        # Plot Segmented Data
        for segment in segments:
            START_END = [x_offset + segment['start'], x_offset + segment['end']]
            SEG_CN = [segment['segmented'], segment['segmented']]
            # plot segs
            plt.plot(
                START_END,
                SEG_CN,
                color=cat_colours[segment['category']],
                linewidth=2
            )
    chr_limits=chr_limits[:-1]
    # Customizing the plot
    plt.ylim(-3, 3)
    # plt.axhline(0, color='black', linestyle='--', linewidth=1)
    for lim in chr_limits:
        plt.axvline(x=lim, color='black', linewidth=0.25)
    # # add gois
    # for gene in gois_plotting:
    #     plt.axvline(x=gene[0], color=gene_colour, linewidth=0.5)
    #     plt.text(gene[0]+500000,3,gene[1],rotation=90, color=gene_colour, fontsize=5)
    plt.margins(x=0) 
    # plt.suptitle(title, fontsize=5) 
    # plt.title(subtitle,fontsize=5)
    # plt.title(title,fontsize=5, loc='left', fontweight='bold')
    # plt.title(title,fontsize=5, loc='left')
    # plt.title(subtitle,fontsize=5, loc='right')
    plt.title(f"{title}\n{subtitle}", fontsize=5, loc='left', pad=2)
    # plt.title('Copy Number Profile')
    plt.xlabel('Chromosome',fontsize=5)
    plt.ylabel('Copy number (log2R)',fontsize=5)
    # set thin axis spine linewidth for all sides
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(0.2)
    plt.tick_params(width=0.2)
    desired_labels = {str(i) for i in range(1, 12)}.union({str(i) for i in range(12, 23, 2)})
    # keep all tick positions but only show text for the desired labels;
    # use empty string for ticks we want to keep but not label
    display_labels = [lab if lab in desired_labels else '' for lab in chr_label]
    ax.set_xticks(midpoints)
    ax.set_xticklabels(display_labels, fontsize=5)
    # plt.xticks(midpoints, chr_label,fontsize=5)
    plt.yticks(fontsize=5)
    # plt.legend(loc='upper right', markerscale=4)
    plt.grid(False)
    plt.tight_layout()  
    # plt.rcParams.update({'font.size': 5})  
    # plt.rc('font', **font)
    # Save plot
    plt.savefig(fig_name, dpi=300)
    plt.close()
        

# calculate chromosome offset for plotting
def calc_chrom_offset(seg_data):
    # estimate chr length
    chrom_ends = {}
    for i in range(len(seg_data)):
        chr = (seg_data[i]['chrom'])
        if chr not in chrom_ends:
            chrom_ends[chr] = [seg_data[i]['end']]
        elif chr in chrom_ends:
            chrom_ends[chr].append(seg_data[i]['end'])
    # estimate offsets
    chromosome_data = {}
    offset = 0
    for chrom in list(dict.fromkeys(chrom_ends)):
        chr_length = max(chrom_ends[chrom])
        start = offset + 1
        end = offset + chr_length
        mid = offset + round(chr_length / 2)
        chromosome_data[chrom] = {
            'length': chr_length,
            'offset': offset,
            'start': start,
            'end': end,
            'mid': mid
        }
        offset += chr_length
    return chromosome_data


def plot_copy_number(absolute_cn_path, log2r_cn_path, cn_fit_path, sample, no_plot_points, outdir):
# def plot_copy_number(absolute_cn_path, log2r_cn_path, cn_fit_path, goi_annots_path, sample, outdir):
    # Color mapping for copy number categories
    cat_colours = {
        'del': '#0072B2',
        'loss': '#56B4E9',
        'neut': '#A9A9A9',
        'gain': '#E69F00',
        'amp': '#D55E00',
        'other': '#D3D3D3'
    }
    # gene_colour = "#009E73"
    # segments
    seg_info = {}
    with open(absolute_cn_path, "r") as file:
        header=next(file)
        for line in file:
            fields = line.strip().split("\t")
            seg_id, chr, start, end, category = fields[3],fields[0],int(fields[1]),int(fields[2]),fields[-1]
            seg_info[seg_id] = {'chrom':chr,
                                'start':start,
                                'end':end,
                                'category':category}
    # bins
    bin_data = []
    segs = {}
    with open(log2r_cn_path, "r") as file:
        header=next(file)
        for line in file:
            fields = line.strip().split("\t")
            chr, position, bin_cn, seg_cn, seg_id = fields[1], (int(fields[2])+int(fields[3])-1)/2, float(fields[-3]), float(fields[-1]), fields[-2]
            cn_cat = seg_info[seg_id]['category']
            if seg_id not in segs:
                segs[seg_id] = seg_cn
            bin_data.append({'chrom':chr,
                             'pos':float(position),
                             'binned': bin_cn,
                             'category': cn_cat})
    # downsample bin data for plotting
    if len(bin_data) > no_plot_points:
        bin_data = random.sample(bin_data, no_plot_points)
    # add relative segment copy number to segment data
    seg_data = []
    for seg_id in list(dict.fromkeys(segs)):
        seg_info[seg_id]['segmented'] = segs[seg_id]
        seg_data.append(seg_info[seg_id])
    # read in sample fitting metrices
    with open(cn_fit_path, "r") as file:
        header=next(file)
        for line in file:
            fields = line.strip().split("\t")
            meta = {'sample': sample,
                    'purity': fields[4],
                    'fitted_purity':fields[0],
                    'ploidy': fields[1],
                    'gof': fields[2],
                    'coverage': fields[7],
                    'pag': fields[8]}
    # # prep genes of interest
    # if goi_annots_path is not None:
    #     goi_list = []
    #     with open(goi_annots_path, "r") as file:
    #         header=next(file)
    #         for line in file:
    #             fields = line.strip().split("\t")
    #             if fields[-1] != 'neut':
    #                 goi_list.append([fields[0],fields[1],fields[2],fields[3]])
    # estimate chromosome offsets
    offsets = calc_chrom_offset(seg_data)
    # plot data 
    # plotting(bin_data, seg_data, offsets, goi_list, cat_colours, gene_colour, meta, outdir)
    plotting(bin_data, seg_data, offsets, cat_colours, meta, outdir)


### focal region plotting ###
def focal_plotting(outdir, meta, offsets, cat_colours, region_data, roi_coords, density_points):
    # Prepare figure
    sample,roi_name = meta['sample'],meta['roi_name']
    title = f"{sample} | region: {roi_name}"
    subtitle = f"CN deviation={meta['cnsep']}; T-statistic={meta['T']}; $\it{{P}}${meta['pval']}"
    fig_name = f"{outdir}/{sample}_focal_analysis_{roi_name}.pdf"
    fig = plt.figure(figsize=(3.5, 1.25))
    gs = fig.add_gridspec(1, 2, width_ratios=[6, 1], wspace=0)
    p_main = fig.add_subplot(gs[0])
    p_side = fig.add_subplot(gs[1], sharey=p_main)
    # plot roi cross
    if roi_coords != None:
        p_main.axvline(x=roi_coords[0]+offsets[roi_coords[2]]['offset'], color=cat_colours['roi'], linewidth=0.75)
        p_main.axhline(y=roi_coords[1], color=cat_colours['roi'], linewidth=0.75)
        if density_points != None:
            p_side.axhline(y=roi_coords[1], color=cat_colours['roi'], linewidth=0.75)
    # prepare chromosome input fpr plotting
    chromosomes = list(dict.fromkeys(offsets))
    chr_label = [x.replace('chr','') for x in chromosomes]
    midpoints = []
    chr_limits = []
    # Plot each chromosome
    for chrom in chromosomes:
        # offsets
        x_offset = offsets[chrom]['offset']
        midpoints.append(offsets[chrom]['mid'])
        chr_limits.append(offsets[chrom]['end'])            
        # Filter data for the current chromosome
        bins = [d for d in region_data if d['chrom'] == chrom]
        # Plot Region log2R Data
        for cat, col in cat_colours.items():
            subset = [d for d in bins if d['category'] == cat]
            POS = [x_offset + d['pos'] for d in subset]
            CN = [d['binned'] for d in subset]
            # scatter bins
            p_main.scatter(
                POS,
                CN,
                color=col,
                s=1 if cat == 'background' else 4,
                alpha=0.75
            )
    chr_limits=chr_limits[:-1]
    # Customizing the plot
    minval,maxval = min([d['binned'] for d in region_data])-1, max([d['binned'] for d in region_data])+1
    plt.ylim(minval,maxval)
    # plt.ylim((maxval*-1),maxval)
    for lim in chr_limits:
        p_main.axvline(x=lim, color='black', linewidth=0.25)
    p_main.set_xlabel('Chromosome',fontsize=5)
    p_main.set_ylabel('Copy number (log2R)',fontsize=5)
    # p_main.set_xticks(midpoints, chr_label,fontsize=5)
    # only label chromosomes 1-10 and even chromosomes >=12 (12,14,...,22)
    
    # desired_labels = [str(i) for i in range(1, 12)] + [str(i) for i in range(13, 23, 2)]
    # selected = [(pos, lab) for pos, lab in zip(midpoints, chr_label) if lab in desired_labels]
    # if selected:
    #     sel_pos, sel_lab = zip(*selected)
    #     p_main.set_xticks(sel_pos)
    #     p_main.set_xticklabels(sel_lab, fontsize=5)
    # else:
    #     p_main.set_xticks(midpoints)
    #     p_main.set_xticklabels(chr_label, fontsize=5)
    
    desired_labels = {str(i) for i in range(1, 12)}.union({str(i) for i in range(12, 23, 2)})
    # keep all tick positions but only show text for the desired labels;
    # use empty string for ticks we want to keep but not label
    display_labels = [lab if lab in desired_labels else '' for lab in chr_label]
    p_main.set_xticks(midpoints)
    p_main.set_xticklabels(display_labels, fontsize=5)

    # p_main.text(0.0, 1.05, title, fontsize=5, ha='left', va='bottom', transform=p_main.transAxes)
    # p_main.text(1.0, 1.05, subtitle, fontsize=5, ha='right', va='bottom', transform=p_main.transAxes)
    p_main.set_title(f"{title}\n{subtitle}", fontsize=5, loc='left', pad=2)
    p_main.tick_params(axis='y', labelsize=5)
    p_main.grid(False)
    p_main.margins(y=0)
    p_main.margins(x=0.01)
    p_main.set_xlim(p_main.get_xlim()[0], p_main.get_xlim()[1])
    # change all spines to thinner lines
    for axis in ['top','bottom','left','right']:
        p_main.spines[axis].set_linewidth(0.2)
    p_main.tick_params(width=0.2)
    ## add in side density plot ##
    if density_points != None:
        y_vals, density = zip(*density_points)
        p_side.fill_betweenx(y_vals, 0, density, color=cat_colours['background'], alpha=0.6)
        p_side.plot(density, y_vals, color=cat_colours['background'], linewidth=0.5)
        p_side.set_xlabel('Density', fontsize=5)
    p_side.tick_params(axis='y', left=False, labelleft=False)
    p_side.tick_params(axis='x', bottom=False, labelbottom=False)
    p_side.spines['right'].set_visible(False)
    p_side.spines['top'].set_visible(False)
    p_side.spines['bottom'].set_visible(False)
    p_side.spines['left'].set_linewidth(0.2)
    # plot asthetics
    plt.yticks(fontsize=5)
    plt.xticks(fontsize=5)
    plt.margins(x=0) 
    plt.tight_layout()  
    # plt.title(title,fontsize=5, loc='center')
    # plt.title(subtitle,fontsize=5)
    # Save plot
    plt.savefig(fig_name, dpi=300)
    plt.close()

def label_pvalue(pval):
    pval = float(pval)
    if pval <= 0.0001:
        # plab = '≤0.0001****'
        plab = '≤0.0001'
    elif pval <= 0.001:
        # plab = '≤0.001***'
        plab = '≤0.001'
    elif pval <= 0.01:
        # plab = '≤0.01**'
        plab = '≤0.01'
    elif pval <= 0.05:
        # plab = '≤0.05*'
        plab = '≤0.05'
    elif pval > 0.05:
        plab = '=NS'
        # plab = '>0.05'
    return plab
    
def plot_focal_results(focal_cn, focal_out, sample, outdir):
    # prepare input data
    ## meta data
    meta = {'sample': sample,
            'roi_name': focal_out[3],
            'T': round(float(focal_out[8]),3) if focal_out[8] not in ['None', 'NA'] else focal_out[8],
            'pval': label_pvalue(float(focal_out[9])) if focal_out[9] not in ['None', 'NA'] else focal_out[9],
            'cnsep': round(float(focal_out[10]),3) if focal_out[10] not in ['None', 'NA'] else focal_out[10]}
    ## cn data
    roi_coords = ( (int(focal_out[1])+int(focal_out[2])-1)/2 , float(focal_out[4]), focal_out[0]) if focal_out[4] != 'NA' else None
    region_data = []
    for line in focal_cn:
        chr, start, end, position, bin_cn = line[1], int(line[2]), int(line[3]), (int(line[2])+int(line[3])-1)/2, float(line[-1])
        label = line[-2] if line[-2] != "background" else None
        bin_cat = line[-2] if line[-2] == "background" else "roi"
        region_data.append({'chrom':chr,
                            'pos':float(position),
                            'start':start,
                            'end':end,
                            'binned': bin_cn,
                            'label': label,
                            'category': bin_cat})
    # chromosome offsets
    offsets = calc_chrom_offset(region_data)
    # colours
    cat_colours = {
        'background': '#A9A9A9',
        'roi': '#CC79A7'
    }
    # get density data for side plot
    bg_cn = [float(x[-1]) for x in focal_cn if x[-2] == "background"]
    try:
        dens_x,dens_y = cnfitter.r_density_default(bg_cn, n=512)
        density_points = [(x, y) for x, y in zip(dens_x,dens_y)]
    except:
        print(f"density could not be computed. values of n = {len(bg_cn)} sampled.")
        density_points = None
    # plot
    focal_plotting(outdir, meta, offsets, cat_colours, region_data, roi_coords, density_points)
    
