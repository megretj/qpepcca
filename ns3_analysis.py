
folder_path = "./ns3_results/19.07.23/hybla/"
# Get a list of all files in the folder
file_list = os.listdir(folder_path)
average_500 = {}
average_100 = {}
average_30 = {}
# Loop through each file and save the average
for i,filename in enumerate(file_list):
    # Search for the pattern in the filename
    match = re.search(r"0\.(\d+?)-", filename)
    if not match: 
        continue
    error = float("0."+ match.group(1))
    rtt = re.search(r"-(\d+?)\.",filename).group(1)
    data = np.loadtxt(f'./ns3_results/19.07.23/hybla/{filename}')
    if rtt == "250":
        average_500[error] = np.average(data[:,1])/0.5
    if rtt == "50":
        average_100[error] = np.average(data[:,1])/0.1
    if rtt == "15":
        average_30[error] = np.average(data[:,1])/0.03
average_500=dict(sorted(average_500.items()))
average_100=dict(sorted(average_100.items()))
average_30=dict(sorted(average_30.items()))