import os
import subprocess
import random
import json
from Bio.Seq import Seq
from wtforms import validators
from flask_wtf import FlaskForm, CSRFProtect
from werkzeug.utils import secure_filename
from wtforms import StringField, BooleanField, SubmitField, SelectField
from wtforms.validators import DataRequired, Length
from flask import Flask, render_template, request, redirect, url_for, flash, abort, current_app, jsonify
from flask_wtf.file import FileField
from wtforms import StringField, TextAreaField, SubmitField
import re
from tool.window_folding_based_selection import get_potential_windows_scores
from tool.switch_generator import SwitchGenerator
from google.cloud import storage
from tool.server_utils import process_file_stream
import os

# Initialize the Flask app
app = Flask(__name__, template_folder='tool/templates')
csrf = CSRFProtect(app)
csrf.init_app(app)
app.config['SECRET_KEY'] = os.urandom(24)
app.config['MAX_CONTENT_LENGTH'] = 210 * 1024 * 1024  # 210 MB
app.config['UPLOAD_EXTENSIONS'] = ['.fasta']

client = storage.Client()
bucket_name = 'protech_bucket'
bucket = client.get_bucket(bucket_name)

class InputForm(FlaskForm):
    email = StringField("Email", [DataRequired()], render_kw={"id": "email"})
    target_seq = StringField("RNA sequence", validators=[DataRequired()], render_kw={"id": "gene"})
    user_trigger_bool = BooleanField("Got a known trigger?", render_kw={"id": "user_trigger_bool"})
    trigger = StringField("Input Trigger (23 nucleotides)", render_kw={"id": "trigger"}, validators=[validators.Length(min=23, max=23, message="Trigger must be 23 nucleotides long")])
    reporter_gene = StringField("Reporter Gene", render_kw={"id": "reporter_gene"}, validators=[DataRequired()])
    cell_type = SelectField("Organism Type",
                            choices=[('Prokaryote', 'Prokaryote'), ('Eukaryote', 'Eukaryote'),
                                     ('Homo sapiens', 'Homo sapiens')], render_kw={"id": "cell_type"})
    file = FileField('Transcripts File', render_kw={"id": "file"})
    submit = SubmitField("Submit")

# Home page route
@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')


@app.route('/form', methods=['GET', 'POST'])
def user_data_getter():
    email = None
    target_seq = None
    user_trigger_bool = False
    trigger = None
    reporter_gene = None
    cell_type = None
    file = None

    input_form = InputForm()
    if input_form.validate_on_submit():
        try:
            # Get the data from the form
            email = input_form.email.data
            target_seq = input_form.target_seq.data.upper()
            trigger = input_form.trigger.data.upper()
            reporter_gene = input_form.reporter_gene.data.upper()
            cell_type = input_form.cell_type.data
            user_trigger_bool = input_form.user_trigger_bool.data
            uploaded_file = input_form.file.data


            # Process the file
            s_email = str(email) if email else "EMPTY"
            s_target_seq = str(target_seq) if target_seq else "EMPTY"
            s_trigger = str(trigger) if trigger else "EMPTY"
            s_reporter_gene = str(reporter_gene) if reporter_gene else "EMPTY"
            s_cell_type = str(cell_type) if cell_type else "EMPTY"
            s_user_trigger_bool = 'True' if user_trigger_bool else 'EMPTY'

            if uploaded_file:
                try:
                    s_file_dict = process_file_stream(uploaded_file)
                except Exception as e:
                    flash(f"Error processing file: {e} check validate sequences")
                    s_file_dict = None
                s_file_dict = json.dumps(s_file_dict) if s_file_dict else "EMPTY"
            else:
                s_file_dict = "EMPTY"

            # TODO: ADD SUBPROCESS ID
            subprocess.run(['python', 'tool/generate.py', s_email, s_target_seq, s_trigger, s_reporter_gene, s_cell_type,
                            s_user_trigger_bool, s_file_dict],
                           text=True)

            print("Subprocess ran successfully")
            input_form.email.data = ''
            input_form.target_seq.data = ''
            input_form.trigger.data = ''
            input_form.reporter_gene.data = ''
            input_form.user_trigger_bool.data = ''
            input_form.cell_type.data = ''
            input_form.file.data = ''
            flash('Form submitted successfully. Job accepted.')
        except Exception as e:
            flash(f"Error processing form: {e}")

    return render_template('form.html', input_form=input_form)

# Error Handling
@app.errorhandler(404)
def page_not_found(e):
    return render_template("404.html"), 404

@app.errorhandler(500)
def internal_error(e):
    return render_template("500.html"), 500


if __name__ == '__main__':
    app.run(debug=True)
