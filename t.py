import requests

headers = {'Content-Type': 'application/x-www-form-urlencoded'}

data = 'and(keyword(%221000565%22),or(keyword(%22METUNv1_00054%22),keyword(%22metunv1_00054%22),keyword(%22METUNv1_00055%22),keyword(%22metunv1_00055%22),keyword(%22METUNv1_00056%22),keyword(%22metunv1_00056%22),keyword(%22METUNv1_00057%22),keyword(%22metunv1_00057%22),keyword(%22METUNv1_00058%22),keyword(%22metunv1_00058%22),keyword(%22METUNv1_00059%22),keyword(%22metunv1_00059%22),keyword(%22METUNv1_00060%22),keyword(%22metunv1_00060%22),keyword(%22METUNv1_00061%22),keyword(%22metunv1_00061%22),keyword(%22METUNv1_00062%22),keyword(%22metunv1_00062%22),keyword(%22METUNv1_00063%22),keyword(%22metunv1_00063%22),keyword(%22METUNv1_00064%22),keyword(%22metunv1_00064%22),keyword(%22METUNv1_00065%22),keyword(%22metunv1_00065%22),keyword(%22METUNv1_00066%22),keyword(%22metunv1_00066%22),keyword(%22METUNv1_00067%22),keyword(%22metunv1_00067%22),keyword(%22METUNv1_00068%22),keyword(%22metunv1_00068%22),keyword(%22METUNv1_00069%22),keyword(%22metunv1_00069%22),keyword(%22METUNv1_00070%22),keyword(%22metunv1_00070%22),keyword(%22METUNv1_00071%22),keyword(%22metunv1_00071%22),keyword(%22METUNv1_00072%22),keyword(%22metunv1_00072%22),keyword(%22METUNv1_00073%22),keyword(%22metunv1_00073%22)))&limit(1)'

response = requests.post('https://patricbrc.org/api/genome_feature', headers=headers, data=data)
response.raise_for_status()
print(dir(response))
