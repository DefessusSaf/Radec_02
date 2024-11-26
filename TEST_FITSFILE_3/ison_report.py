import xml.etree.ElementTree as ET
import os


def parse_data_with_filename(raw_data, sensor_id=10121):
    """
    Преобразует данные в структуру, добавляя информацию о файле.
    
    :param raw_data: Список строк данных (имя файла, временная метка, координаты).
    :param sensor_id: Идентификатор источника данных.
    :return: Список словарей, соответствующих объектам.
    """
    # Извлекаем имя файла и временную метку
    file_name = raw_data[0].split(": ")[1].strip()  # Первой строкой - имя файла
    utc = raw_data[1].strip()  # Вторая строка - временная метка
    
    # Остальные строки - координаты объектов
    coordinates = raw_data[2:]
    parsed_data = []
    
    # Преобразуем данные координат в структуру
    for idx, coord in enumerate(coordinates, start=1):
        ra, dec = coord.strip().split()  # Разделяем RA и Dec
        # ra, dec = format_conservation(ra, dec)
        parsed_data.append({
            "file": file_name,
            "sensor": sensor_id,
            "id": idx,  # Генерация ID объекта
            "utc": utc,
            "ra_j2000": ra,
            "dec_j2000": dec,
            "mag": None,
            "suspicious": False
        })
    
    return parsed_data


def create_ison_report(data, output_file):
    """
    Создает отчет в формате XML из переданных данных.
    
    :param data: Список данных, которые нужно записать.
    :param output_file: Путь к файлу для сохранения отчета.
    """
    root = ET.Element('data')
    
    for entry in data:
        meas = ET.SubElement(root, 'meas')
        ET.SubElement(meas, 'file').text = entry['file']
        ET.SubElement(meas, 'sensor').text = str(entry['sensor'])
        ET.SubElement(meas, 'id').text = str(entry['id'])
        ET.SubElement(meas, 'utc').text = entry['utc']
        ET.SubElement(meas, 'ra_j2000').text = entry['ra_j2000']
        ET.SubElement(meas, 'dec_j2000').text = entry['dec_j2000']
        if entry['mag'] is not None:
            ET.SubElement(meas, 'mag').text = str(entry['mag'])
        ET.SubElement(meas, 'suspicious').text = str(entry['suspicious']).lower()
    
    tree = ET.ElementTree(root)
    tree.write(output_file, encoding='utf-8', xml_declaration=True)


# Пример исходных данных
# raw_data = [
#     "#File: IRS1B_Rocket_00002.fits",
#     "2024-10-31T12:50:57.4433554",
#     "0:55:43.15 71:38:03.01",
#     "0:56:42.90 71:16:22.13"
# ]

def process_file_in_dir(dir):
    all_data = []
    for filename in os.listdir(dir):
        file_path = os.path.join(dir, filename)
        if os.path.isfile(file_path):
            with open(file_path, "r") as file:
                raw_data = file.readlines()
            raw_data = [line.strip() for line in raw_data if line.strip() and not line.startswith("#")]
            
            parsed_data = parse_data_with_filename(raw_data)
            all_data.extend(parsed_data)
    
    return all_data

dir = "PROCESS_FILE"

# Обрабатываем данные
all_data  = process_file_in_dir(dir)
# Создаем отчет в XML формате
create_ison_report(all_data, "ison_report.xml")

print("Отчет успешно создан!")