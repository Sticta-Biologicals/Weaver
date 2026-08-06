import datetime
from django import forms
from .models import Plasmid
from .models import Primer
from .models import GlycerolStock
from .custom.primer_access import visible_primers_for_user
from .custom.standards import assembly_standards
from organization.views import get_projects_where_member_can
from organization.views import get_projects_where_member_can_any


class PlasmidNameInput(forms.Form):
    plasmid_name = forms.CharField(max_length=50, min_length=1)


class GlycerolQRInput(forms.Form):
    glycerol_qr_id = forms.CharField()


class BlastSequenceInput(forms.Form):
    def __init__(self, project_choices, *args, **kwargs):
        super(BlastSequenceInput, self).__init__(*args, **kwargs)
        self.fields['project'].choices = project_choices

    project = forms.ChoiceField(label="Project to search in", choices=())
    fasta_sequence = forms.CharField(label="Fasta Text Input (Preferred)", widget=forms.Textarea(attrs={}), required=False)
    fasta_file = forms.FileField(label="Fasta File Input", required=False)
    short_blast = forms.BooleanField(label="Use short input BLAST parameters?", required=False)


class L0SequenceInput(forms.Form):
    l0_sequence_input = forms.CharField(widget=forms.Textarea, label="Sequence input")
    # Todo append None
    assembly_standard_slug = forms.ChoiceField(choices=tuple([(index, assembly_standard['name']) for index, assembly_standard in assembly_standards.items()]),
                                required=True, label="Ligation standard")
    l0_oh_5 = forms.ChoiceField(choices=tuple([(oh_slug, oh['name'] + " [" + oh['oh'] + "]") for index, assembly_standard in assembly_standards.items() for oh_slug, oh in assembly_standard['ohs']['l0'].items()]),
                                required=True, label="L0 OH 5'")
    l0_oh_3 = forms.ChoiceField(choices=tuple([(oh_slug, oh['name'] + " [" + oh['oh'] + "]") for index, assembly_standard in assembly_standards.items() for oh_slug, oh in assembly_standard['ohs']['l0'].items()]),
                                required=True, label="L0 OH 3'")
    enzyme = forms.CharField(required=True, widget=forms.HiddenInput())


class FastaFileInput(forms.ClearableFileInput):
    allow_multiple_selected = True


class FastaFileField(forms.FileField):
    def clean(self, data, initial=None):
        single_file_clean = super().clean
        if isinstance(data, (list, tuple)):
            return [single_file_clean(d, initial) for d in data]
        return single_file_clean(data, initial)


class FastaAlignForm(forms.Form):
    def __init__(self, *args, **kwargs):
        super(FastaAlignForm, self).__init__(*args, **kwargs)
        self.fields['fasta_sequence'].widget.attrs.update({'class': 'form-control', 'rows': 10, 'placeholder': '>query\nACGT...'})
        self.fields['fasta_file'].widget.attrs.update({'class': 'form-control', 'accept': '.fa,.fas,.fasta,.txt'})
        self.fields['alignment_view_mode'].widget.attrs.update({'class': 'form-select'})
        self.fields['save_clustal_file'].widget.attrs.update({'class': 'form-check-input'})

    fasta_sequence = forms.CharField(widget=forms.Textarea(attrs={}), required=False)
    alignment_view_mode = forms.ChoiceField(
        label="View mode",
        choices=(("combined", "Together"), ("individual", "One at a time")),
        initial="combined",
    )
    save_clustal_file = forms.BooleanField(required=False)
    fasta_file = FastaFileField(label="Fasta File", required=False, widget=FastaFileInput(attrs={"multiple": True}))


class MultipleFileInput(forms.ClearableFileInput):
    allow_multiple_selected = True


class MultipleFileField(forms.FileField):
    def clean(self, data, initial=None):
        single_file_clean = super().clean
        if isinstance(data, (list, tuple)):
            return [single_file_clean(d, initial) for d in data]
        return single_file_clean(data, initial)


class SangerAlignForm(forms.Form):
    save_clustal_file = forms.BooleanField(required=False)
    label = forms.CharField(label="Run label", required=False, max_length=200)
    notes = forms.CharField(label="Notes", required=False, widget=forms.Textarea(attrs={"rows": 2}))
    sanger_files = MultipleFileField(label="Sanger files", required=False, widget=MultipleFileInput(attrs={
        "multiple": True,
        "accept": ".ab1,.phd.1,.seq,.fa,.fas,.fasta",
        "id": "id_sanger_files",
    }))
    ab1 = forms.FileField(label="AB1 File", required=False)

    def clean(self):
        cleaned_data = super().clean()
        if not self.files.getlist("sanger_files") and not self.files.get("ab1"):
            raise forms.ValidationError("Upload at least one .ab1, .phd.1, or .seq file.")
        return cleaned_data


class DateInput(forms.DateInput):
    input_type = 'date'


class GstockCreateForm(forms.ModelForm):
    class Meta:
        model = GlycerolStock
        fields = '__all__'
        exclude = ('project',)
        widgets = {
            'created_on': DateInput()
        }


class GstockEditForm(forms.ModelForm):
    class Meta:
        model = GlycerolStock
        fields = '__all__'
        widgets = {
            'created_on': DateInput()
        }


class PlasmidCreateForm(forms.ModelForm):
    def __init__(self, *args, **kwargs):
        member = kwargs.pop('user')
        super(PlasmidCreateForm, self).__init__(*args, **kwargs)
        user_visible_plasmids = Plasmid.objects.filter(project__in=get_projects_where_member_can_any(member)).order_by('name')
        # self.fields['backbone'].queryset = user_visible_plasmids.filter(type=1)
        # self.fields['inserts'].queryset = user_visible_plasmids.filter(type=0)

    class Meta:
        model = Plasmid
        fields = ['name', 'selectable_markers', 'sequence', 'backbone', 'inserts', 'intended_use', 'type', 'level',
                  'description', 'project', 'ligation_state', 'created_on']


class PlasmidEditForm(forms.ModelForm):
    def __init__(self, *args, **kwargs):
        member = kwargs.pop('user')
        super(PlasmidEditForm, self).__init__(*args, **kwargs)
        user_visible_plasmids = Plasmid.objects.filter(project__in=get_projects_where_member_can_any(member)).order_by('name')
        self.fields['project'].queryset = get_projects_where_member_can(member, ['a', 'w'])
        # self.fields['backbone'].queryset = user_visible_plasmids.filter(type=1)
        # self.fields['inserts'].queryset = user_visible_plasmids.filter(type=0)

    class Meta:
        model = Plasmid
        fields = ['name', 'selectable_markers', 'sequence', 'backbone', 'inserts', 'intended_use', 'type', 'level',
                  'description', 'created_on', 'project', 'ligation_state', 'reference_sequence']


class PlasmidValidationForm(forms.ModelForm):
    class Meta:
        model = Plasmid
        fields = ['ligation_state', 'working_colony',
                  'colonypcr_state', 'colonypcr_date', 'colonypcr_observations',
                  'digestion_state', 'digestion_date', 'digestion_observations',
                  'sequencing_state', 'sequencing_date', 'sequencing_observations',
                  'sequencing_observations', 'sequencing_clustal_file']
        widgets = {
            'colonypcr_date': DateInput(),
            'digestion_date': DateInput(),
            'sequencing_date': DateInput()
        }


class DigestForm(forms.Form):
    enzymes = forms.CharField(widget=forms.HiddenInput(attrs={'id': 'digest_enzymes'}))


class PCRForm(forms.Form):
    def __init__(self, *args, **kwargs):
        user = kwargs.pop('user', None)
        super(PCRForm, self).__init__(*args, **kwargs)
        if user is not None:
            queryset = visible_primers_for_user(user).order_by('name')
            self.fields['primer_f'].queryset = queryset
            self.fields['primer_r'].queryset = queryset

    primer_f = forms.ModelChoiceField(queryset=Primer.objects.all(), to_field_name="id", label="Primer F",
                                      required=False)
    primer_r = forms.ModelChoiceField(queryset=Primer.objects.all(), to_field_name="id", label="Primer R",
                                      required=False)
    primer_f_seq = forms.CharField(label="Primer F sequence", required=False)
    primer_r_seq = forms.CharField(label="Primer R sequence", required=False)


class PrimerBatchUploadForm(forms.Form):
    def __init__(self, *args, **kwargs):
        super(PrimerBatchUploadForm, self).__init__(*args, **kwargs)
        self.fields['fasta_file'].widget.attrs.update({'class': 'form-control', 'accept': '.fa,.fas,.fasta,.txt'})
        self.fields['name_source'].widget.attrs.update({'class': 'form-select'})
        self.fields['default_direction'].widget.attrs.update({'class': 'form-select'})
        self.fields['update_existing'].widget.attrs.update({'class': 'form-check-input'})

    fasta_file = forms.FileField(label="FASTA file")
    update_existing = forms.BooleanField(label="Update existing primers with the same name", required=False)
    name_source = forms.ChoiceField(
        label="Primer name source",
        choices=(("id", "FASTA ID"), ("description", "Full FASTA description")),
        initial="id",
    )
    default_direction = forms.ChoiceField(
        label="Direction when F/R is not in the name",
        choices=(("f", "Forward"), ("r", "Reverse"), ("", "Leave empty")),
        initial="f",
        required=False,
    )


class ServicesPCRForm(forms.Form):
    def __init__(self, *args, **kwargs):
        user = kwargs.pop('user', None)
        super(ServicesPCRForm, self).__init__(*args, **kwargs)
        if user is not None:
            primers = visible_primers_for_user(user).order_by('name')
            self.fields['primer_f'].queryset = primers.filter(fwd_or_rev='f')
            self.fields['primer_r'].queryset = primers.filter(fwd_or_rev='r')
        for field in self.fields.values():
            field.widget.attrs.update({'class': 'form-control'})
        self.fields['primer_f'].widget.attrs.update({'class': 'form-select'})
        self.fields['primer_r'].widget.attrs.update({'class': 'form-select'})

    primer_f = forms.ModelChoiceField(queryset=Primer.objects.none(), to_field_name="id", label="Forward primer")
    primer_r = forms.ModelChoiceField(queryset=Primer.objects.none(), to_field_name="id", label="Reverse primer")
    min_product_size = forms.IntegerField(label="Minimum product size", min_value=1, initial=100)
    max_product_size = forms.IntegerField(label="Maximum product size", min_value=1, required=False)

    def clean(self):
        cleaned_data = super().clean()
        min_product_size = cleaned_data.get('min_product_size')
        max_product_size = cleaned_data.get('max_product_size')
        if min_product_size and max_product_size and max_product_size < min_product_size:
            raise forms.ValidationError("Maximum product size must be greater than or equal to minimum product size.")
        return cleaned_data


class BatchPrintsForm(forms.Form):
    LABEL_TYPES = (
        ("plasmids", "Plasmids"),
        ("glycerolstocks", "Glycerol stocks"),
    )

    def __init__(self, *args, **kwargs):
        super(BatchPrintsForm, self).__init__(*args, **kwargs)
        self.fields["label_type"].widget.attrs.update({"class": "form-select"})
        self.fields["identifiers"].widget.attrs.update({
            "class": "form-control",
            "rows": 8,
            "placeholder": "One ID, code, QR, or name per line",
        })
        self.fields["date"].widget.attrs.update({"class": "form-control"})
        self.fields["concentration"].widget.attrs.update({"class": "form-control"})

    label_type = forms.ChoiceField(label="Label type", choices=LABEL_TYPES)
    identifiers = forms.CharField(label="Items", widget=forms.Textarea)
    date = forms.DateField(widget=DateInput(), initial=datetime.date.today)
    concentration = forms.FloatField(label="Concentration (ng/ul)", min_value=0.01, initial=1)


class MsaUploadAb1FilesForm(forms.Form):
    ab1_file_1 = forms.FileField(label="AB1 File 1")
    ab1_file_2 = forms.FileField(label="AB1 File 2", required=False)


class MsaChromatosStep2Form(forms.Form):
    from1 = forms.IntegerField(label="From Chromato 1")
    to1 = forms.IntegerField(label="To Chromato 1")
    from2 = forms.IntegerField(label="From Chromato 2", required=False)
    to2 = forms.IntegerField(label="To Chromato 2", required=False)
    sequence1 = forms.CharField(widget=forms.HiddenInput)
    sequence2 = forms.CharField(widget=forms.HiddenInput, required=False)
    target = forms.ModelChoiceField(queryset=Plasmid.objects.all(), to_field_name="id")


class MsaUploadFastaFileForm(forms.Form):
    fasta_text = forms.CharField(widget=forms.Textarea)
    fasta_file = forms.FileField(label="Fasta file")


class PlasmidLabel(forms.Form):
    date = forms.DateField(widget=DateInput(), initial=datetime.date.today)
    colony = forms.CharField()
    concentration = forms.CharField()
