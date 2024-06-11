import os
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import simple_norm

if __name__ == '__main__':
    # Load the data
    image_path = './output/end_to_end_demo_VIS.fits'
    if not os.path.exists(image_path):
        raise FileNotFoundError(f'File {image_path} not found')
    image = fits.getdata(image_path)

    # Plot the data
    plt.figure()
    norm = simple_norm(image, 'asinh', asinh_a=0.1, vmin=-0.01, vmax=0.8)
    plt.imshow(image, norm=norm, origin='lower', cmap='Greys_r')
    plt.title('Euclid VIS image')
    plt.colorbar()
    plt.show()
