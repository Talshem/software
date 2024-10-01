from .server_utils import process_file_stream
from .utils.seqeunce_consts import M_CHERRY_ORF, GFP_GENE
from .server_utils import validate_sequence_bool


import os
import subprocess
import json
from flask_wtf import FlaskForm, CSRFProtect
from wtforms.validators import DataRequired
from flask import Flask, render_template, flash, request
from flask_wtf.file import FileField
from wtforms import StringField, SubmitField, BooleanField, SelectField, ValidationError, validators
from werkzeug.utils import secure_filename
#from google.cloud import storage
from flask import session

# Initialize the Flask app
app = Flask(__name__, template_folder='templates')
csrf = CSRFProtect(app)
csrf.init_app(app)

app.config['SECRET_KEY'] = 'shaniqua'
app.config['MAX_CONTENT_LENGTH'] = 210 * 1024 * 1024  # 210 MB
app.config['UPLOAD_EXTENSIONS'] = ['.fasta']
"""
client = storage.Client()
bucket_name = 'protech_bucket'
bucket = client.get_bucket(bucket_name)
"""
def validate_trigger_length(form, field):
    organism_type = form.cell_type.data
    if not form.user_trigger_bool.data:
        return True

    if not field.data:
        raise ValidationError("Trigger is required.")
    if not validate_sequence_bool(field.data):
        raise ValidationError("Invalid sequence.")

    trigger_len = len(field.data)
    trigger_bool = form.user_trigger_bool.data
    if trigger_bool:
        if organism_type == 'E.coli' and trigger_len != 23:
            raise ValidationError("Trigger must be 23 nucleotides long for prokaryotes.")
        elif organism_type in ['Saccharomyces cerevisiae', 'Homo sapiens'] and trigger_len != 30:
            raise ValidationError("Trigger must be 30 nucleotides long for eukaryotes or humans.")
    return True

def validate_sequence(form, field):
    if form.user_trigger_bool.data:
        return True
    if not field.data:
        raise ValidationError("Sequence is required.")
    if not validate_sequence_bool(field.data):
        raise ValidationError("Invalid sequence.")
    return True

def validate_reporter(form, field):
    if not field.data:
        if form.optional_reporter.data == 'None':
            raise ValidationError("Sequence is required or choose from optional reporters.")
        else:
            return True
    else:
        if form.optional_reporter.data != 'None':
            raise ValidationError("Sequence is required or choose from optional reporters, can't input both.")
        else:
            return True
    if not validate_sequence_bool(field.data):
        raise ValidationError("Invalid sequence.")

    return True


def validate_file_format(form, filed):
    if filed.data:
        format = f'.{str(secure_filename(filed.data.filename)).split(".")[-1]}'
        if format in ['.fasta', '.fa', '.fastq', '.gb', '.gbk', '.embl', '.phy', '.phylip', '.aln', '.nex',
                      '.stockholm', '.tab',
                      '.qual', '.abi', '.sff', '.xml', '.ig']:
            return True
        else:
            raise ValidationError("Invalid format.")
    else:
        return True


class InputForm(FlaskForm):
    email = StringField("Email", [DataRequired()], render_kw={"id": "email"})
    target_seq = StringField("RNA sequence", render_kw={"id": "gene"}, validators=[validate_sequence])
    cell_type = SelectField("Organism", choices=[('E.coli', 'E.coli'), ("Saccharomyces cerevisiae", 'Saccharomyces cerevisiae'),
                                     ('Homo sapiens', 'Homo sapiens')], render_kw={"id": "cell_type"})
    optional_reporter = SelectField("Optional Reporters",
                            choices=[('None', 'None'), ('GFP', 'GFP'),('mCherry','mCherry')], render_kw={"id": "optional_reporter"})
    user_trigger_bool = BooleanField("Got a known trigger?", render_kw={"id": "user_trigger_bool"})
    trigger = StringField("Known Trigger", render_kw={"id": "trigger"}, validators=[validate_trigger_length])
    reporter_gene = StringField("Reporter Gene Sequence", render_kw={"id": "reporter_gene"}, validators=[validate_reporter])
    file = FileField('Transcripts File', render_kw={"id": "file"}, validators=[validate_file_format])
    submit = SubmitField("Submit")


# Home page route
@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@app.route('/success', methods=['GET'])
def success_page():
    return render_template('success.html')


@app.route('/form', methods=['GET', 'POST'])
def user_data_getter():
    email = None
    target_seq = None
    user_trigger_bool = False
    trigger = None
    reporter_gene = None
    cell_type = None
    file = None
    reporter_option = None

    input_form = InputForm()
    if input_form.validate_on_submit():
        try:
            email = input_form.email.data
            if trigger:
                trigger = input_form.trigger.data.upper()
            else:
                target_seq = input_form.target_seq.data.upper()

            reporter_gene = input_form.reporter_gene.data.upper()
            reporter_option = input_form.optional_reporter.data
            if reporter_option != 'None':
                reporter_gene = reporter_option
            else:
                if reporter_gene == 'GFP':
                    reporter_gene = GFP_GENE
                else:
                    reporter_gene = M_CHERRY_ORF

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
                format = f'.{str(secure_filename(uploaded_file.filename)).split(".")[-1]}'
                try:
                    s_file_dict = process_file_stream(uploaded_file, format)
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
