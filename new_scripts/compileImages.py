#
# Note, this is a python3 script
#
# Usage: python3 compileImages.py <pattern> <conditions.txt> <Ntiles> <Image Directory>


import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def autoContrast(im):
    '''Port of ImageJ's autocontrast function. Returns a clipped image
       with min and max intensities set by the algorithm
       Adapted from Kota Miura, <http://wiki.cmci.info/documents/120206pyip_cooking/python_imagej_cookbook>'''

    auto_threshold = 5000

    pixel_count = len(im)**2

    limit = int(pixel_count/10)
    histogram = np.unique(im,return_counts=True)

    threshold = int(pixel_count/auto_threshold)

    count = 0
    for i in np.arange(len(histogram[0])):
        count = histogram[1][i]
        if count > limit:
            count = 0
        if count > threshold:
            break

    vmin = histogram[0][i]

    for i in np.arange(len(histogram[0])-1,0,-1):
        count = histogram[1][i]
        if count > limit:
            count = 0
        if count> threshold:
            break
    vmax = histogram[0][i]

    clipped_im = np.clip(im,vmin,vmax)

    return clipped_im


### BEGIN SCRIPT ###

if __name__ == '__main__':


    ## Check input arguments
    if len(sys.argv) == 3:
        print('Using default <Ntiles> = 18 and <Image Directory> = Images/')
        pattern = sys.argv[1]
        conditions_file = sys.argv[2]
        Ntiles = 18
        out_dir = 'Images'

    elif len(sys.argv) == 5:
        pattern = sys.argv[1]
        conditions_file = sys.argv[2]
        Ntiles = int(sys.argv[3])
        out_dir = sys.argv[4]
    else:
        print('Usage: python compileImages.py <pattern> <conditions.txt> <Ntiles> <Image Directory>')
        sys.exit(0)


    ## Deal with input paramters
    conditions = []
    with open(conditions_file, 'r') as f:
        for line in f:
            condition = line.strip()
            conditions += [condition]

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)


    ## Get a list of all images
    image_dirs = [dir for dir in os.listdir('.') if pattern in dir]

    images = []
    for image_dir in image_dirs:
        try:
            images += [*[image_dir+'/'+f for f in os.listdir(image_dir)]]
        except NotADirectoryError:
            print(f"NotADirectoryError, {image_dir}")
            continue



    ## Loop through tiles and make plots
    for tile_no in np.arange(Ntiles)+1:
        green_tiles = np.sort([tile for tile in images if re.search(f'tile{tile_no}[^\d]+green',tile)])
        red_tiles = np.sort([tile for tile in images if re.search(f'tile{tile_no}[^\d]+red',tile)])

        Ncond = max(len(green_tiles),len(red_tiles))

        print(Ncond)
        if Ncond == 0:
            print('skip')
            continue
        else:
            print(Ncond)

        fig = plt.figure(figsize=(8,11))
        subfigs = gridspec.GridSpec(1, 2, figure=fig)

        title = fig.suptitle(f'Images for Tile {tile_no}',fontsize=24)

        red_axs = subfigs[0].subgridspec(int(np.ceil(Ncond/2)),2)
        green_axs = subfigs[1].subgridspec(int(np.ceil(Ncond/2)),2)

        for i,c in zip((0,1),('Red','Green')):
            subtitle = fig.add_subplot(subfigs[0,i])
            subtitle.set_title(f'{c} Channel Images\n',fontsize=20)
            subtitle.set_frame_on(False)
            subtitle.axis('off')

        for i,condition in zip(range(Ncond),conditions):
            try:
                red_im = autoContrast(plt.imread(red_tiles[i]))
                red_ax = fig.add_subplot(red_axs[np.unravel_index(i,(int(np.ceil(Ncond/2)),2),order='F')])
                red_ax.imshow(red_im,cmap='gist_gray')

                red_ax.set_title(condition,fontsize=16)
                red_ax.set_xticklabels('')
                red_ax.set_yticklabels('')
                red_ax.set_xticks([])
                red_ax.set_yticks([])

            except:
                print(f'no red image for condition {condition}')

            try:
                green_im = autoContrast(plt.imread(green_tiles[i]))
                green_ax = fig.add_subplot(green_axs[np.unravel_index(i,(int(np.ceil(Ncond/2)),2),order='F')])
                green_ax.imshow(green_im,cmap='gist_gray')

                green_ax.set_title(condition,fontsize=16)
                green_ax.set_xticklabels('')
                green_ax.set_yticklabels('')
                green_ax.set_xticks([])
                green_ax.set_yticks([])
            except:
                print(f'no red image for condition {condition}')


        plt.tight_layout()
        fig.savefig(f'{out_dir}/Tile_{tile_no}.png',facecolor='white',bbox_extra_artists=[title])
        fig.clf()
        plt.close()
