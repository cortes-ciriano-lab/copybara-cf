import matplotlib.pyplot as plt


# plotting
def plotting(bin_data, seg_data, offsets, cat_colours, meta):
    # Prepare figure
    title = f"{meta['sample']}"
    subtitle = f"PAG={meta['pag']}; purity={meta['purity']}; fitted purity={meta['fitted_purity']}; fitted ploidy={meta['ploidy']}; goodness of fit={meta['gof']}; coverage={meta['coverage']}x"
    fig_name = f"{title}_copy_number_plot.pdf"
    plt.figure(figsize=(7, 2))
    # prepare chromosome input fpr plotting
    chromosomes = list(dict.fromkeys(offsets))
    chr_label = [x.replace('chr','') for x in chromosomes]
    midpoints = []
    # chr_limits = [offsets['chr1']['start']]
    chr_limits = []
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
        plt.axvline(x=lim, color='black', linewidth=0.5)
    plt.margins(x=0) 
    # plt.suptitle(title, fontsize=7) 
    # plt.title(subtitle,fontsize=6)
    # plt.title(title,fontsize=7, loc='left', fontweight='bold')
    plt.title(title,fontsize=7, loc='left')
    plt.title(subtitle,fontsize=7, loc='right')
    # plt.title('Copy Number Profile')
    plt.xlabel('Chromosome',fontsize=6)
    plt.ylabel('Copy number (log2R)',fontsize=6)
    plt.xticks(midpoints, chr_label,fontsize=5)
    # plt.legend(loc='upper right', markerscale=4)
    plt.grid(False)
    plt.tight_layout()  
    # plt.rcParams.update({'font.size': 6})  
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


def plot_copy_number(absolute_cn_path, log2r_cn_path, cn_fit_path, sample):
    # Color mapping for copy number categories
    cat_colours = {
        'del': '#0072B2',
        'loss': '#56B4E9',
        'neut': '#A9A9A9',
        'gain': '#E69F00',
        'amp': '#D55E00',
        'other': '#D3D3D3'
    }
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
    # estimate chromosome offsets
    offsets = calc_chrom_offset(seg_data)
    # plot data 
    plotting(bin_data, seg_data, offsets, cat_colours, meta)
