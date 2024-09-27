# Email configuration
import os
import smtplib
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from io import BytesIO

import pandas as pd

SMTP_SERVER = 'smtp.gmail.com'
SMTP_PORT = 587
EMAIL_ADDRESS = 'erlichnet57@gmail.com'
EMAIL_PASSWORD = "bzkb uqbp ykgs dkqc"

def send_email_with_attachment(to_email, results_df):
    subject = 'ProTech WebTool Report'
    body = 'Please find the attached Excel file with the results of your analysis.'
    msg = MIMEMultipart()
    msg['From'] = EMAIL_ADDRESS
    msg['To'] = to_email
    msg['Subject'] = subject
    msg.attach(MIMEText(body, 'plain'))

    # Convert DataFrame to Excel in-memory using BytesIO
    excel_buffer = BytesIO()
    results_df.to_excel(excel_buffer, index=False, engine='xlsxwriter')  # Write the DataFrame to buffer
    excel_buffer.seek(0)  # Move the cursor back to the beginning of the buffer


    # Attach the Excel file in memory to the email
    part = MIMEBase('application', 'octet-stream')
    part.set_payload(excel_buffer.read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', 'attachment', filename='results.xlsx')
    msg.attach(part)

    # Send the email
    server = smtplib.SMTP(SMTP_SERVER, SMTP_PORT)
    server.starttls()  # Secure the connection
    server.login(EMAIL_ADDRESS, EMAIL_PASSWORD)  # Login to the email account
    server.sendmail(EMAIL_ADDRESS, to_email, msg.as_string())  # Send the email
    server.quit()

    return

"""
def send_email(results_df):
    # Connect to the server
    server = smtplib.SMTP(SMTP_SERVER, SMTP_PORT)
    server.starttls()  # Secure the connection
    server.login(EMAIL_ADDRESS, EMAIL_PASSWORD)  # Login to the email account

    # Send the email
    server.send_message(results_df)
    server.quit()
    print("Email sent successfully!")
"""

# Example usage
if __name__ == "__main__":

    file = '/Users/netanelerlich/PycharmProjects/webTool/utils/path_to_your_excel_file.xlsx'
    email = 'erlichnet57@gmail.com'
    df = pd.read_csv(file)
    df  = pd.DataFrame(df)
    send_email_with_attachment(email, df)
