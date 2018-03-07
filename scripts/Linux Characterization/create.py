import csv

csv_file = "stacks.csv"
in_txt = csv.reader(open('stacks.txt', "rb"), delimiter = '#')
out_csv = csv.writer(open(csv_file, 'wb'))
out_csv.writerows(in_txt)

csv_file = "conf_pucker.csv"
i_txt = csv.reader(open('conf_pucker.txt', "rb"), delimiter = '#')
out_csv = csv.writer(open(csv_file, 'wb'))
out_csv.writerows(i_txt)


csv_file = "hbonds.csv"
txt = csv.reader(open('hbonds.txt', "rb"), delimiter = '#')
out_csv = csv.writer(open(csv_file, 'wb'))
out_csv.writerows(txt)

csv_file = "bp.csv"
txt = csv.reader(open('bp.txt', "rb"), delimiter = '#')
out_csv = csv.writer(open(csv_file, 'wb'))
out_csv.writerows(txt)

csv_file = "detailedbp.csv"
txt = csv.reader(open('detailedbp.txt', "rb"), delimiter = '#')
out_csv = csv.writer(open(csv_file, 'wb'))
out_csv.writerows(txt)