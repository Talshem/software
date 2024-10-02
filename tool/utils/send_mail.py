# Email configuration
import smtplib
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from io import BytesIO
import base64



SMTP_SERVER = 'smtp.gmail.com'
SMTP_PORT = 587
SENDER_EMAIL_ADDRESS = 'igemtau2024@gmail.com'
EMAIL_PASSWORD = "eefg btov iwsp xzgt"


def send_email_with_attachment(results_df, to_email):
    subject = 'ToREC WebTool Analysis Report'
    body = """
    Dear User,

    We are pleased to inform you that your analysis has been completed successfully. 

    Please find attached the detailed report in Excel format, which contains the results of your toehold switch design and analysis based on free energy calculations and competition modeling. We hope this data provides valuable insights for your project.

    If you have any questions or need further assistance, feel free to reach out.

    Thank you for using the WebTool, and we look forward to assisting you in future analyses.

    Best regards,
    The ProTech Team
    """

    msg = MIMEMultipart()
    msg['From'] = SENDER_EMAIL_ADDRESS
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
    server.login(SENDER_EMAIL_ADDRESS, EMAIL_PASSWORD)  # Login to the email account
    server.sendmail(SENDER_EMAIL_ADDRESS, to_email, msg.as_string())  # Send the email
    server.quit()

    return


# Example usage
if __name__ == "__main__":
    pass


