import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from astropy.coordinates import Angle
import astropy.units as u
import os
import glob

def is_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sqrt(np.sum((points - median) ** 2, axis=-1))
    med_abs_deviation = np.median(diff)
    modified_z_score = 0.6745 * diff / med_abs_deviation
    return np.where(modified_z_score > thresh)[0]


def load_data(file_path):
    data =  np.genfromtxt(file_path, unpack=True)
    return data


def choose_fits_file():
    log_file = 'TMP/processing_log.txt'
    if not os.path.exists(log_file):
        raise FileNotFoundError("Log file does not exist.")

    with open(log_file, 'r') as f:
        last_file = f.readline().strip()

    if not last_file:
        raise ValueError("No file name found in the log file.")

    last_file_base = os.path.basename(last_file)
    # Извлекаем суффикс (все символы до расширения)
    last_file_suffix = last_file_base.split('.')[0]
    print(f"Extracted suffix: {last_file_suffix}")

    fits_files = glob.glob('TMP/*.fits') + glob.glob('TMP/*.fit')
    
    # Печать всех найденных файлов
    print("Found FITS files:")
    for f in fits_files:
        print(f)

    matching_file = None
    for fits_file in fits_files:
        file_name = os.path.basename(fits_file)
        # Печать для отладки
        print(f"Checking file: {file_name}")
        # Проверка наличия суффикса в имени файла
        if last_file_suffix in file_name:
            matching_file = fits_file
            break

    if not matching_file:
        raise FileNotFoundError(f"No matching FITS file found for suffix {last_file_suffix}")

    print(f"Selected FITS file: {matching_file}")

    return matching_file, last_file_base


def preprocess_data(X, Y, ERRX, ERRY, A, B, XMIN, YMIN, XMAX, YMAX, TH, FLAG, FLUX, x_min, x_max, y_min, y_max):
    y_mask = (Y >= y_min) & (Y <= y_max)
    x_mask = (X >= x_min) & (X <= x_max)
    mask = x_mask & y_mask
    return (X[mask], Y[mask], ERRX[mask], ERRY[mask], A[mask], B[mask], 
            XMIN[mask], YMIN[mask], XMAX[mask], YMAX[mask], TH[mask], FLAG[mask], FLUX[mask])


def compute_elongation(A, B):
    return A / B


def compute_errores(ERRX, ERRY):
    erroreX = np.sqrt(1/ERRX)
    erroreY = np.sqrt(1/ERRY)
    # print(f"erroreX: {erroreX} erroreY: {erroreY}")
    return erroreX, erroreY


def convert(ra_deg, dec_deg):
    ra_angle = Angle(ra_deg, unit=u.deg)
    ra_hms = ra_angle.to_string(unit=u.hour, sep=':', precision=2)
    dec_angle = Angle(dec_deg, unit=u.deg)
    dec_dms = dec_angle.to_string(unit=u.deg, sep=':', precision=2)
    return ra_hms, dec_dms


def save_results(coords, fits_filename, base_filename, X, Y, ERRX, ERRY, A, B, XMIN, YMIN, XMAX, YMAX):
    output_dir = 'PROCESS_FILE'
    os.makedirs(output_dir, exist_ok=True)
    
    txt_filename = os.path.join(output_dir, f'{base_filename}.txt')
    with open(txt_filename, 'w') as f:
        f.write(f"File: {base_filename}\n")
        
        with fits.open(fits_filename) as hdul:
            header = hdul[0].header
            date_obs = header.get('DATE-OBS', '00000')
            f.write(f"{date_obs}\n")
        with fits.open(fits_filename) as hdul:
            wcs = WCS(hdul[0].header)
        sky_coords = wcs.pixel_to_world(X, Y)
        for coord, x, y, errx, erry, a, b, xmin, ymin, xmax, ymax in zip(
                sky_coords, X, Y, ERRX, ERRY, A, B, XMIN, YMIN, XMAX, YMAX):
            ra_deg = coord.ra.deg
            dec_deg = coord.dec.deg
            ra_hms, dec_dms = convert(ra_deg, dec_deg)
            f.write(f"{ra_hms} {dec_dms} {x} {y} {errx} {erry} {a} {b} {xmin} {ymin} {xmax} {ymax}\n")


def star_observation(X, Y, ERRX, ERRY, A, B, XMIN, YMIN, XMAX, YMAX, fits_filename, base_filename):
    ELONG = compute_elongation(A, B)
    # print(f"ELONG values: {ELONG}")
    
    elong_mask = ELONG > 2.2
    # print(f"ELONG mask: {elong_mask}")
    
    X_filtered = X[elong_mask]
    Y_filtered = Y[elong_mask]
    if len(X_filtered) == 0:
        print("No data passed the ELONG filter.")
        return

    erroreX, erroreY = compute_errores(ERRX, ERRY)
    print(f"Filtered X: {X_filtered}, Y: {Y_filtered}")

    with fits.open(fits_filename) as hdul:
        wcs = WCS(hdul[0].header)

    try:
        sky_coords = wcs.pixel_to_world(X_filtered, Y_filtered)
    except Exception as e:
        print(f"WCS conversion error: {e}")
        return

    coords_filtered = []
    for coord in sky_coords:
        ra_deg = coord.ra.deg
        dec_deg = coord.dec.deg
        ra_hms, dec_dms = convert(ra_deg, dec_deg)
        coords_filtered.append((ra_hms, dec_dms))
        print(f"RA: {ra_hms}, DEC: {dec_deg}")

    save_results(coords_filtered, fits_filename, base_filename, 
                 X_filtered, Y_filtered, erroreX[elong_mask], erroreY[elong_mask], 
                 A[elong_mask], B[elong_mask], XMIN[elong_mask], YMIN[elong_mask], 
                 XMAX[elong_mask], YMAX[elong_mask])

    return coords_filtered

def main():
    DIR = 'TMP/'
    fn = 'k1-impTEST.fts.sx'
    try:
        fits_filename, base_filename = choose_fits_file()
    except (FileNotFoundError, ValueError) as e:
        print(e)
        return

    X, Y, ERRX, ERRY, A, B, XMIN, YMIN, XMAX, YMAX, TH, FLAG, FLUX = load_data(os.path.join(DIR, fn))
    # print(f"X: {X}, Y: {Y}, ERRX: {ERRX}, ERRY: {ERRY}, A: {A}, B: {B}, XMIN: {XMIN}, YMIN: {YMIN}, XMAX: {XMAX}, YMAX: {YMAX}, TH: {TH}, FLAG: {FLAG}, FLUX: {FLUX}")
    
    X, Y, ERRX, ERRY, A, B, XMIN, YMIN, XMAX, YMAX, TH, FLAG, FLUX = preprocess_data(X, Y, ERRX, ERRY, A, B, XMIN, YMIN, XMAX, YMAX, TH, FLAG, FLUX, x_min=25, x_max=4785, y_min=15, y_max=3175)
    # print(f"X: {X}, Y: {Y}, ERRX: {ERRX}, ERRY: {ERRY}, A: {A}, B: {B}, XMIN: {XMIN}, YMIN: {YMIN}, XMAX: {XMAX}, YMAX: {YMAX}, TH: {TH}, FLAG: {FLAG}, FLUX: {FLUX}")
    
    star_observation(X, Y, ERRX, ERRY, A, B, XMIN, YMIN, XMAX, YMAX, fits_filename, base_filename)

if __name__ == "__main__":
    main()