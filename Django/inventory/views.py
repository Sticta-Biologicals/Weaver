from django.shortcuts import render
from django.views.decorators.csrf import csrf_exempt
from django.utils.decorators import method_decorator
from django.utils.text import get_valid_filename
from django.shortcuts import redirect

from .models import Plasmid
from .models import GlycerolStock
from .models import RestrictionEnzyme
from .models import Primer
from .models import Box
from .models import Location
from .models import SangerRead
from .models import SangerReadFile
from .models import SangerVariant
from .models import SangerVerificationRun
from .custom.general import CHECK_STATES
from .custom.general import LIGATION_STATES
from .custom.pcr import suggest_pcr_primers
from .custom.pcr import matching_primer_annotations
from .custom.pcr import matching_amplicon_annotations
from .custom.pcr import primer_pair_amplicons
from .custom.pcr import primer_pair_complementarity
from .custom.pcr import select_non_overlapping_amplicons
from .custom.pcr import amplicon_segments
from .custom.pcr import display_primer_name
from .custom.pcr import display_primer_id
from .custom.restriction_digest import DigestConstraints
from .custom.restriction_digest import DEFAULT_MAX_ENZYMES
from .custom.restriction_digest import DEFAULT_MAX_FRAGMENTS
from .custom.restriction_digest import DEFAULT_MIN_BAND_DIFFERENCE_BP
from .custom.restriction_digest import DEFAULT_MIN_BUFFER_ACTIVITY_PERCENT
from .custom.restriction_digest import DEFAULT_MIN_FRAGMENT_SIZE_BP
from .custom.restriction_digest import DEFAULT_MIN_FRAGMENTS
from .custom.restriction_digest import DEFAULT_RESULT_LIMIT
from .custom.restriction_digest import enzymes_with_effective_cuts
from .custom.restriction_digest import normalize_regions
from .custom.restriction_digest import serialize_digest_response
from .custom.primer_access import visible_primers_for_user
from .custom.primer_import import PrimerImportError
from .custom.primer_import import import_primers_from_fasta
from .custom.sanger import alignment_tracks_for_ove
from .custom.sanger import align_read
from .custom.sanger import classify_run
from .custom.sanger import clustal_content
from .custom.sanger import combined_metrics
from .custom.sanger import display_trim_range
from .custom.sanger import process_sanger_files
from .custom.sanger import parse_ab1
from .custom.sanger import read_is_usable
from .custom.sanger import read_metrics_tsv
from .custom.sanger import SangerProcessingParameters
from .custom.sanger import trim_by_quality
from .custom.sanger import variants_csv
from .models import Stats
from .models import PlasmidType
from .models import TableFilter

from organization.decorators import require_current_project_set
from organization.decorators import require_member_can_any_current_project
from organization.decorators import require_member_can_write_or_admin_current_project
from organization.decorators import require_member_can_read_project_of_plasmid
from organization.decorators import require_member_can_read_project_of_gs
from organization.decorators import require_member_can_read_project_of_primer
from organization.decorators import require_member_can_write_or_admin_project_of_plasmid
from organization.decorators import require_member_can_write_or_admin_project_of_gs
from organization.decorators import require_member_can_write_or_admin_project_of_primer
from organization.views import get_current_project_id
from organization.views import get_current_project
from organization.views import on_current_project_member_can_write_or_admin
from organization.views import get_projects_where_member_can_any
from organization.views import member_can_write_or_admin_plasmid
from organization.views import member_can_write_or_admin_gs
from organization.views import get_show_from_all_projects
from organization.views import member_can_write_or_admin_primer
from organization.views import on_project_member_can_any

from .forms import PlasmidValidationForm
from .forms import PlasmidCreateForm
from .forms import PlasmidEditForm
from .forms import PlasmidLabel

from django.http import HttpResponseRedirect

import html
from django.conf import settings
from django.http import HttpResponse, Http404
from django.core.exceptions import ObjectDoesNotExist
from .custom.box import BOX_ROWS
from .custom.box import BOX_COLUMNS
from django.views.generic.edit import UpdateView
from django.views.generic.edit import CreateView
from django.views.generic.edit import DeleteView
from django.urls import reverse
from django.core.exceptions import PermissionDenied
from django.db import transaction
from django.db.models import Q

from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq
from Bio.Seq import reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Restriction import RestrictionBatch
from Bio.Restriction import AllEnzymes
from Bio.Restriction.Restriction_Dictionary import rest_dict
from pyblast import BioBlast
from pyblast.utils import make_linear, make_circular

DEFAULT_AMPLICON_REGION_FLANK_BP = 30

from .forms import DigestForm
from .forms import PCRForm
from .forms import PrimerBatchUploadForm
from .forms import ServicesPCRForm
from .forms import BatchPrintsForm
import json
from io import StringIO
from .forms import SangerAlignForm
from .forms import FastaAlignForm
from .forms import L0SequenceInput
from .forms import BlastSequenceInput
from .forms import GstockCreateForm
from .forms import GstockEditForm
from .forms import GlycerolQRInput
from .forms import PlasmidNameInput

import os
import tempfile
import re
import uuid
import builtins
import Bio
from Bio import AlignIO
import django
from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove
from django.core.files.base import ContentFile
from datetime import datetime
from datetime import date
from django.utils import timezone
from contextlib import contextmanager
from django.http import JsonResponse
from io import StringIO
import requests
from bs4 import BeautifulSoup

import plotly.express as px
import pandas as pd


@contextmanager
def pyblast_open_compat():
    original_open = builtins.open

    def open_compat(file, mode='r', *args, **kwargs):
        if mode == 'rU':
            mode = 'r'
        return original_open(file, mode, *args, **kwargs)

    builtins.open = open_compat
    try:
        yield
    finally:
        builtins.open = original_open


def run_pyblast_compat(callback):
    with pyblast_open_compat():
        return callback()


def get_table_filters(level_from_table_filters, level_to_table_filters):
    pt_table_filters = []
    for pt in PlasmidType.objects.all():
        pt_table_filters.append((pt.name, 't' + str(pt.id), 'success'))

    level_table_filters = []
    for level in range(level_from_table_filters, level_to_table_filters + 1):
        level_table_filters.append(('L' + str(level), 'l' + str(level), 'warning'))

    sw_table_filters = []
    for tf in TableFilter.objects.all():
        filters = []
        color = 'info'
        if tf.color:
            color = tf.color
        for the_filter in tf.options.split(","):
            name, search = the_filter.split("|")
            filters.append((name, search, color))
        sw_table_filters.append(['startswith', tf.name, filters])

    table_filters = [
                        ['all', 'All', [
                            ('All', 'all', 'primary'),
                        ]],
                        ['type', 'Type', pt_table_filters],
                        ['level', 'Level', level_table_filters]
                    ] + sw_table_filters
    return table_filters


def json_serial(obj):
    """JSON serializer for objects not serializable by default json code"""
    if isinstance(obj, (datetime, date)):
        return obj.isoformat()
    raise TypeError("Type %s not serializable" % type(obj))


def restrictionenzyme(request, restrictionenzyme_id):
    try:
        restrictionenzyme_to_detail = RestrictionEnzyme.objects.get(id=restrictionenzyme_id)
    except ObjectDoesNotExist:
        raise Http404
    context = {
        'restrictionenzyme': restrictionenzyme_to_detail,
    }
    return render(request, 'inventory/restrictionenzyme.html', context)


def restrictionenzymes(request):
    context = {
        'restrictionenzymes': RestrictionEnzyme.objects.all(),
    }
    return render(request, 'inventory/restrictionenzymes.html', context)


@require_current_project_set
def glycerolstocks(request):
    show_from_all_projects = get_show_from_all_projects(request)
    if show_from_all_projects:
        glycerolstocks = GlycerolStock.objects.filter(project_id__in=get_projects_where_member_can_any(request.user))
    else:
        glycerolstocks = GlycerolStock.objects.filter(project_id=get_current_project_id(request))

    level_from_table_filters = 0
    level_to_table_filters = 0
    hasGlycerolStocks = False
    for glycerolstock in glycerolstocks:
        hasGlycerolStocks = True
        if glycerolstock.plasmid:
            if glycerolstock.plasmid.level:
                if glycerolstock.plasmid.level > level_to_table_filters:
                    level_to_table_filters = glycerolstock.plasmid.level
                if glycerolstock.plasmid.level < level_from_table_filters:
                    level_from_table_filters = glycerolstock.plasmid.level
    context = {
        'on_current_project_member_can_write_or_admin': on_current_project_member_can_write_or_admin(request),
        'has_glycerolstocks': hasGlycerolStocks,
        'table_filters': get_table_filters(level_from_table_filters, level_to_table_filters),
        'show_from_all_projects': show_from_all_projects
    }
    return render(request, 'inventory/glycerolstocks.html', context)


@require_member_can_read_project_of_gs
def glycerolstock(request, glycerolstock_id):
    try:
        glycerolstock_to_detail = GlycerolStock.objects.get(id=glycerolstock_id)
    except ObjectDoesNotExist:
        raise Http404

    resistantes_human_context = "None"
    if glycerolstock_to_detail.plasmid:
        resistantes_human_context = resistantes_human(glycerolstock_to_detail.plasmid.selectable_markers)

    glycerolstock_to_detail.resistantes_human = resistantes_human_context
    glycerolstock_to_detail.resistantes_strain_human = resistantes_human(
        glycerolstock_to_detail.strain.selectable_markers, True)

    context = {
        'glycerolstock': glycerolstock_to_detail,
        'user_can_edit_gs': member_can_write_or_admin_gs(glycerolstock_to_detail, request.user),
    }
    return render(request, 'inventory/glycerolstock.html', context)


def glycerolstock_qr(request):
    if request.method == 'POST' and 'glycerol_qr_id' in request.POST:
        form = GlycerolQRInput(request.POST)
        if form.is_valid():
            return glycerolstock_from_qr(request, form.cleaned_data['glycerol_qr_id'])
    else:
        context = {
            'glycerol_qr_nput_form': GlycerolQRInput()
        }
    return render(request, 'inventory/glycerolstock_qr.html', context)


def glycerolstock_from_qr(request, glycerolstock_qr_id):
    try:
        glycerolstock_to_detail = GlycerolStock.objects.filter(qr_id=glycerolstock_qr_id)[0]
    except IndexError:
        raise Http404
    except ObjectDoesNotExist:
        raise Http404
    return redirect('glycerolstock', glycerolstock_id=glycerolstock_to_detail.id)


def gstock_check_pos(the_class, the_self, the_form):
    try:
        obj_curr_pos = GlycerolStock.objects.get(box_row=the_form.cleaned_data['box_row'],
                                                 box_column=the_form.cleaned_data['box_column'],
                                                 box=the_form.cleaned_data['box'])
        if obj_curr_pos:
            if obj_curr_pos != the_form.instance:
                the_form.add_error(None, 'Box position not available')
                return super(the_class, the_self).form_invalid(the_form)
        return super(the_class, the_self).form_valid(the_form)
    except GlycerolStock.DoesNotExist:
        return super(the_class, the_self).form_valid(the_form)


def build_box(request, mode, box):
    box_output = {
        'name': box.name,
        'id': box.id
    }
    if mode == 'p':
        glycerolstocks = box.glycerolstock_set.all()
    else:
        if get_show_from_all_projects(request):
            glycerolstocks = box.glycerolstock_set.filter(
                project_id__in=get_projects_where_member_can_any(request.user))
        else:
            glycerolstocks = box.glycerolstock_set.filter(project_id=get_current_project_id(request))

    anyGs = False
    for glycerolstock in glycerolstocks:
        array_pos = str(glycerolstock.box_row) + str(glycerolstock.box_column)
        if not array_pos in box_output:
            box_output[array_pos] = []
        box_output[array_pos].append(glycerolstock)
        anyGs = True

    return anyGs, box_output


def build_boxes(request, mode):
    output = {
        'BOX_ROWS': BOX_ROWS,
        'BOX_COLUMNS': BOX_COLUMNS,
        'locations': []
    }
    for location in Location.objects.all():
        boxes = []
        for box in Box.objects.filter(location=location).order_by('name'):
            anyGs, box_output = build_box(request, mode, box)

            if anyGs or mode == 'p':
                boxes.append(box_output)
        if boxes or mode == 'p':
            output['locations'].append({
                'name': location.name,
                'boxes': boxes
            })
    return output


class GstockEdit(UpdateView):
    model = GlycerolStock
    template_name_suffix = '_update_form'
    form_class = GstockEditForm

    @method_decorator(require_member_can_write_or_admin_project_of_gs)
    def dispatch(self, *args, **kwargs):
        self.extra_context = {
            'user_can_edit_gs': True
        }
        return super().dispatch(*args, **kwargs)

    def form_valid(self, form):
        return gstock_check_pos(GstockEdit, self, form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["collection"] = build_boxes(self.request, 'p')
        context["render_mod"] = 'p'
        return context

    def get_success_url(self, **kwargs):
        return reverse('glycerolstock', args=(self.object.id,)) + '?form_result_glycerolstock_edit_success=true'


class GstockCreate(CreateView):
    model = GlycerolStock
    template_name_suffix = '_create_form'
    form_class = GstockCreateForm

    @method_decorator(require_member_can_write_or_admin_current_project)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def form_valid(self, form):
        form.instance.project = get_current_project(self.request)
        return gstock_check_pos(GstockCreate, self, form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["collection"] = build_boxes(self.request, 'p')
        context["render_mod"] = 'p'
        return context

    def get_success_url(self, **kwargs):
        return reverse('glycerolstock_label', args=(self.object.id,)) + '?form_result_glycerolstock_create_success=true'


class GstockCreatePlasmidDefined(CreateView):
    model = GlycerolStock
    form_class = GstockCreateForm
    template_name_suffix = '_create_form'

    @method_decorator(require_member_can_write_or_admin_current_project)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def form_valid(self, form):
        form.instance.project = get_current_project(self.request)
        return gstock_check_pos(GstockCreatePlasmidDefined, self, form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["collection"] = build_boxes(self.request, 'p')
        context["plasmid_id"] = self.kwargs['pid']
        context["render_mod"] = 'p'
        return context

    def get_success_url(self, **kwargs):
        return reverse('glycerolstock', args=(self.object.id,)) + '?form_result_glycerolstock_create_success=true'


class GstockDelete(DeleteView):
    model = GlycerolStock

    @method_decorator(require_member_can_write_or_admin_current_project)
    def dispatch(self, *args, **kwargs):
        self.extra_context = {
            'user_can_edit_gs': True
        }
        return super().dispatch(*args, **kwargs)

    def get_success_url(self, **kwargs):
        return reverse('glycerolstock_deleted')


def glycerolstock_deleted(request):
    return render(request, 'inventory/glycerolstock_deleted.html')


@require_member_can_read_project_of_gs
def glycerolstock_label(request, glycerolstock_id):
    try:
        glycerolstock_to_label = GlycerolStock.objects.get(id=glycerolstock_id)
    except ObjectDoesNotExist:
        raise Http404

    resistantes_human_context = "None"
    if glycerolstock_to_label.plasmid:
        resistantes_human_context = resistantes_human(glycerolstock_to_label.plasmid.selectable_markers, True)

    glycerolstock_to_label.resistantes_human = resistantes_human_context
    glycerolstock_to_label.resistantes_strain_human = resistantes_human(
        glycerolstock_to_label.strain.selectable_markers, True)

    context = {
        'glycerolstock': glycerolstock_to_label,
        'user_can_edit_gs': member_can_write_or_admin_gs(glycerolstock_to_label, request.user),
    }
    return render(request, 'inventory/glycerolstock_label.html', context)


def glycerolstock_box(request, box_id):
    try:
        the_box = Box.objects.get(id=box_id)
    except ObjectDoesNotExist:
        raise Http404

    anyGs, box = build_box(request, 'n', the_box)
    context = {
        'collection': {
            'BOX_ROWS': BOX_ROWS,
            'BOX_COLUMNS': BOX_COLUMNS,
        },
        'render_mod': 'n',
        'box': box
    }
    return render(request, 'inventory/glycerolstock_box.html', context)


def glycerolstock_boxes(request):
    context = {
        'collection': build_boxes(request, 'n'),
        'render_mod': 'n',
        'show_from_all_projects': get_show_from_all_projects(request)
    }
    return render(request, 'inventory/glycerolstock_boxes.html', context)


@require_current_project_set
def plasmids(request, render_html=None):
    show_from_all_projects = get_show_from_all_projects(request)
    if show_from_all_projects:
        plasmids = Plasmid.objects.filter(project_id__in=get_projects_where_member_can_any(request.user))
    else:
        plasmids = Plasmid.objects.filter(project_id=get_current_project_id(request))
    level_from_table_filters = 0
    level_to_table_filters = 0
    hasPlasmids = False
    for plasmid in plasmids:
        hasPlasmids = True
        plasmid.refc = plasmid.recommended_enzyme_for_create()
        if plasmid.level:
            if plasmid.level > level_to_table_filters:
                level_to_table_filters = plasmid.level
            if plasmid.level < level_from_table_filters:
                level_from_table_filters = plasmid.level
    context = {
        'on_current_project_member_can_write_or_admin': on_current_project_member_can_write_or_admin(request),
        'table_filters': get_table_filters(level_from_table_filters, level_to_table_filters),
        'has_plasmids': hasPlasmids,
        'RESTRICTION_ENZYMES': RestrictionEnzyme.objects.all,
        'show_from_all_projects': show_from_all_projects
    }
    return render(request, 'inventory/plasmids.html', context)


def resistantes_human(selectable_markers, short=False):
    resistantes_human_return = []
    if selectable_markers:
        for resistance in selectable_markers.all():
            if short:
                resistantes_human_return.append(str(resistance.three_letter_code))
            else:
                resistantes_human_return.append(resistance.name + " (" + str(resistance.three_letter_code) + ")")
            continue
        return " / ".join(resistantes_human_return)
    else:
        return "None"


def recursive_plasmid_build(plasmid):
    if plasmid:
        recursive_result = []
        derivatives = plasmid.get_insert_of()
        if derivatives:
            for derivative in derivatives:
                recursive_result.append((derivative,) + plasmid_create_from_inserts(derivative))
                iter_result = recursive_plasmid_build(derivative)
                if iter_result:
                    recursive_result += iter_result
            
            if len(recursive_result):
                return recursive_result



@require_member_can_read_project_of_plasmid
def plasmid(request, plasmid_id):
    try:
        plasmid_to_detail = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    plasmid_to_detail.refc = plasmid_to_detail.recommended_enzyme_for_create()
    plasmid_to_detail.insert_of = plasmid_to_detail.get_insert_of()
    plasmid_to_detail.backbone_of = plasmid_to_detail.get_backbone_of()

    context = {
        'plasmid': plasmid_to_detail,
        'resistantes_human': resistantes_human(plasmid_to_detail.selectable_markers),
        'CHECK_STATES': CHECK_STATES,
        'RESTRICTION_ENZYMES': RestrictionEnzyme.objects.all(),
        'user_can_edit_plasmid': member_can_write_or_admin_plasmid(plasmid_to_detail, request.user)
    }

    if request.method == 'POST' and 'propagate' in request.POST:
        context['plasmid_propagate_results'] = recursive_plasmid_build(plasmid_to_detail)

    if request.method == 'POST' and 'l0_sequence_input' in request.POST:
        form = L0SequenceInput(request.POST)
        if form.is_valid():
            context['plasmid_create_result'] = plasmid_create_from_inserts(plasmid_to_detail, insert_seq=form.cleaned_data['l0_sequence_input'],
                                        the_re_name=request.POST.get('enzyme'))
        else:
            context['plasmid_create_result'] = ("Bad inputs", "danger")

    if request.method == 'POST' and 'create_from_parts' in request.POST and 'enzyme' in request.POST:
        if plasmid_to_detail.level == 0:
            form = L0SequenceInput(initial={'enzyme': request.POST.get('enzyme')})
            return render(request, 'inventory/plasmid.html',
                          {'L0SequenceInputForm': form, 'plasmid': plasmid_to_detail})
        elif plasmid_to_detail.level == -1:
            form = L0SequenceInput(initial={'enzyme': request.POST.get('enzyme')})
            return render(request, 'inventory/plasmid.html',
                          {'L_1SequenceInputForm': form, 'plasmid': plasmid_to_detail})
        else:
            if 'enzyme' in request.POST:
                context['plasmid_create_result'] = plasmid_create_from_inserts(plasmid_to_detail,
                                            the_re_name=request.POST['enzyme'])
            else:
                context['plasmid_create_result'] = ("No enzyme selected", "danger")

    if request.method == 'GET' and 'ac' in request.GET:
        context['plasmid_create_result'] = plasmid_create_from_inserts(plasmid_to_detail)

    if request.method == 'POST' and 'params' in request.POST:
        context['plasmid_create_result'] = ("Plasmid create wizard is complete.", "success")

    if plasmid_to_detail.public_visibility:
        context['public_url'] = request.build_absolute_uri(reverse('plasmid_public', args=(plasmid_to_detail.id, )))

    # in case of update or never computed
    context['plasmid_update_computed_size_result'] = plasmid_update_computed_size(plasmid_to_detail)

    return render(request, 'inventory/plasmid.html', context)


def plasmid_public(request, plasmid_id):
    try:
        plasmid_to_detail = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        return render(request, 'inventory/general_error.html', {'error': 'Plasmid does not exists or is not public available.'})

    if plasmid_to_detail.public_visibility:
        return render(request, 'inventory/general_error.html', {'error': 'Public system is under construction'})
    else:
        return render(request, 'inventory/general_error.html', {'error': 'Plasmid does not exists or is not public available.'})


def plasmid_create_from_inserts(plasmid_to_build, insert_seq=None, the_re_name=None):
    plasmid_record = plasmid_record_from_inserts(plasmid_to_build, insert_seq, the_re_name)
    if plasmid_record[0]:
        plasmid_record_final = plasmid_record[1]
        plasmid_record_final.name = plasmid_record_final.name.replace(" ", "_")
        plasmid_to_build.sequence.save(plasmid_to_build.name + ".gb", ContentFile(plasmid_record[1].format("gb")))
        return "Plasmid sequence built from backbone / insert data", "success", True
    else:
        return plasmid_record[1], "danger", False


def plasmid_record_from_inserts(plasmid_to_build, insert_seq, the_re_name=None):
    inserts = plasmid_to_build.inserts.all()
    if plasmid_to_build.level < 1 and not insert_seq:
        return False, "Level 0 or -1 plasmids required a defined insert"
    if plasmid_to_build.level > 1 and not inserts and not plasmid_to_build.backbone:
        return False, "No backbone or inserts defined"

    if not the_re_name:
        the_re_name = plasmid_to_build.recommended_enzyme_for_create()

    the_re = None
    for enzyme in AllEnzymes:
        if enzyme.__name__.lower() == the_re_name.lower():
            the_re = enzyme

    if not the_re:
        return False, "Unable to find restriction enzyme"

    cut_length = abs(the_re.ovhg)

    final_record = None

    # L0 & L-1 plasmids
    if plasmid_to_build.level == 0 or plasmid_to_build.level == -1:
        if plasmid_to_build.backbone:
            seqio = seqio_get(plasmid_to_build.backbone)
            if seqio[0]:
                backbone_record = seqio[1]
                hits = the_re.search(backbone_record.seq, linear=False)
                if len(hits) != 2:
                    return False, "Unable to find two restriction sites on " + plasmid_to_build + " backbone"
                # check direction
                if backbone_record.seq[hits[0]-2-len(the_re.site):hits[0]-2].upper() == the_re.site:
                    # looking inwards, rotate
                    backbone_record = extract_circular_region(backbone_record, hits[1]-2-len(the_re.site))
                    hits = the_re.search(backbone_record.seq, linear=False)

                rec_insert = SeqRecord(
                    Seq(insert_seq.lower()),
                    id="oh_h",
                    annotations={"molecule_type": "DNA"}
                )
                rec_insert.features.append(SeqFeature(
                    FeatureLocation(0, len(insert_seq)),
                    type="misc_feature",
                    qualifiers={
                        'ApEinfo_label': plasmid_to_build.name,
                        'label': plasmid_to_build.name,
                        'locus_tag': plasmid_to_build.name,
                    }
                ))

                final_record = rec_insert
                final_record += extract_circular_region(backbone_record, hits[0]-1, length=hits[1]-hits[0])

                final_record.id = str(plasmid_to_build.id)
                final_record.name = plasmid_to_build.name
                final_record.description = plasmid_to_build.description
                final_record.annotations = {"molecule_type": "DNA", "topology": "circular"}

                return True, final_record

            return False, "Error reading backbone sequence file"
        return False, "Plasmid backbone is required for L0 & L(-1) plasmids"

    # >= L1 plasmids
    else:
        parts = []

        if plasmid_to_build.backbone:
            seqio = seqio_get(plasmid_to_build.backbone)
            if seqio[0]:
                backbone_record = seqio[1]
                hits = the_re.search(backbone_record.seq, linear=False)
                if len(hits) != 2:
                    return False, "Unable to find two restriction sites on " + plasmid_to_build + " backbone"
                # check direction
                re_sites_positions = re_site_positions(backbone_record.seq, the_re)
                if hits[0] < re_sites_positions[0] and hits[1] > re_sites_positions[1]:
                    # Invert indexes
                    hits = [hits[1], hits[0]]
                parts.append({
                    "name": plasmid_to_build.backbone.name,
                    "oh5": backbone_record.seq[hits[0]-1:hits[0]-1+cut_length],
                    "oh3": backbone_record.seq[hits[1]-1:hits[1]-1+cut_length],
                    "part": extract_circular_region(backbone_record, hits[0], end=hits[1]),
                })

        if inserts:
            for insert in inserts:
                seqio = seqio_get(insert)
                if seqio[0]:
                    insert_record = seqio[1]
                    hits = the_re.search(insert_record.seq, linear=False)
                    if len(hits) != 2:
                        return False, "Unable to find two restriction sites on " + str(insert) + " insert"
                    # check direction
                    re_sites_positions = re_site_positions(insert_record.seq, the_re)
                    if hits[0] < re_sites_positions[0] and hits[1] > re_sites_positions[1]:
                        # Invert indexes
                        hits = [hits[1], hits[0]]

                    part = {
                        "name": str(insert),
                        "oh5": insert_record.seq[hits[0]-1:hits[0]-1+cut_length],
                        "oh3": insert_record.seq[hits[1]-1:hits[1]-1+cut_length],
                        "part": extract_circular_region(insert_record, hits[0], end=hits[1])
                    }
                    parts.append(part)

        # join parts
        if parts:
            joined = []
            parts_to_iterate = parts.copy()
            if not plasmid_to_build.backbone:
                # first one is a insert
                final_record = parts[0]["part"]
                joined.append(parts[0]["name"] + " (" + str(parts[0]["oh5"]) + "-" + str(parts[0]["oh5"]) + ")")
                parts_to_iterate.pop(0)
            final_oh = parts[0]["oh5"]
            next_oh = parts[0]["oh3"]

            while final_oh != next_oh:
                changed = False
                for idx, part in enumerate(parts_to_iterate):
                    if part["oh5"].upper() == next_oh.upper():
                        if final_record:
                            final_record += part["part"]
                        else:
                            final_record = part["part"]
                        next_oh = part["oh3"]
                        changed = True
                        joined.append(part["name"] + " (" + str(part["oh5"]) + "-" + str(part["oh5"]) + ")")
                        parts_to_iterate.pop(idx)
                        break
                if not changed:
                    return False, "Unable to join parts. Joined = " + " & ".join(joined) + "."

            if plasmid_to_build.backbone:
                # add backbone at last
                final_record += parts[0]["part"]
                joined.append(parts[0]["name"] + " (" + str(parts[0]["oh5"]) + "-" + str(parts[0]["oh5"]) + ")")
                parts_to_iterate.pop(0)

            if len(parts_to_iterate):
                return False, "Parts not participating in the assembly: " + ' / '.join([part["name"] for part in parts_to_iterate])

            final_record.id = str(plasmid_to_build.id)
            final_record.name = plasmid_to_build.name
            final_record.description = plasmid_to_build.description
            final_record.annotations = {"molecule_type": "DNA", "topology": "circular"}

            return True, final_record
        return False, "No parts to assemble"


def re_site_positions(seq, enzyme):
    start = 0
    site_positions = []
    seq = seq.upper()
    site = enzyme.site

    while True:
        start = seq.find(site, start)
        if start == -1:
            break
        site_positions.append(start)
        start += len(enzyme.site)

    start = 0
    site = reverse_complement(site)

    while True:
        start = seq.find(site, start)
        if start == -1:
            break
        site_positions.append(start)
        start += len(enzyme.site)

    return sorted(site_positions)


def extract_circular_region(record, start, length=None, end=None):
    if end:
        if end < start:
            end = len(record.seq) + end
    if not length:
        if end:
            length = end-start
        else:
            length = len(record.seq)-start

    range_start = start - 1
    range_end = start + length - 1
    extended_seq = record.seq + record.seq
    circular_subseq = extended_seq[range_start:range_end]

    # Adjust features
    new_features = []
    for feature in record.features:
        old_feature_start = feature.location.start
        old_feature_end = feature.location.end
        if old_feature_end < range_start or range_end < old_feature_start:
            old_feature_start += len(record.seq)
            old_feature_end += len(record.seq)
            if old_feature_end < range_start or range_end < old_feature_start:
                continue
        feature_start = max(old_feature_start, range_start) - range_start
        feature_end = min(old_feature_end, range_end) - range_start
        if feature_start == feature_end:
            continue
        if feature_end > len(record.seq):
            # goes through origin
            new_features.append(SeqFeature(CompoundLocation([
                FeatureLocation(feature_start, len(record.seq)),
                FeatureLocation(0, feature_end - len(record.seq))
            ], 'join'),  type='rep_origin', qualifiers=feature.qualifiers))
        else:
            new_features.append(SeqFeature(FeatureLocation(feature_start, feature_end), type=feature.type, qualifiers=feature.qualifiers))

    record.seq = circular_subseq
    record.features = new_features

    return record


@require_member_can_read_project_of_plasmid
@require_member_can_write_or_admin_current_project
def plasmid_duplicate(request, plasmid_id):
    try:
        plasmid_to_duplicate = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404
    context = {
        'plasmid': plasmid_to_duplicate,
        'user_can_edit_plasmid': member_can_write_or_admin_plasmid(plasmid_to_duplicate, request.user)
    }
    if 'plasmid_name' in request.POST:
        form = PlasmidNameInput(request.POST)
        context['form'] = form
        if form.is_valid():
            new_plasmid = Plasmid(
                name=form.cleaned_data['plasmid_name'],
                created_on=datetime.now().date(),
                colonypcr_state=1,
                digestion_state=1,
                level=plasmid_to_duplicate.level,
                backbone=plasmid_to_duplicate.backbone,
                type=plasmid_to_duplicate.type,
                sequencing_state=0,
                project=plasmid_to_duplicate.project,
            )
            new_plasmid.save()
            new_plasmid.selectable_markers.add(*plasmid_to_duplicate.selectable_markers.all())
            new_plasmid.inserts.add(*plasmid_to_duplicate.inserts.all())
            return redirect('plasmid', plasmid_id=new_plasmid.id)
        return render(request, 'inventory/plasmid_duplicate.html', context)
    else:
        context['form'] = PlasmidNameInput()
        return render(request, 'inventory/plasmid_duplicate.html', context)


@require_member_can_read_project_of_plasmid
def plasmid_from_qr(request, plasmid_id):
    try:
        plasmid_to_detail = Plasmid.objects.filter(qr_id=plasmid_id)[0]
    except ObjectDoesNotExist:
        raise Http404
    context = {
        'plasmid': plasmid_to_detail,
    }
    return render(request, 'inventory/plasmid.html', context)


class PlasmidEdit(UpdateView):
    model = Plasmid
    form_class = PlasmidEditForm
    template_name_suffix = '_update_form'

    @method_decorator(require_member_can_write_or_admin_project_of_plasmid)
    def dispatch(self, *args, **kwargs):
        self.extra_context = {
            'user_can_edit_plasmid': True
        }
        return super().dispatch(*args, **kwargs)

    def get_form_kwargs(self):
        kwargs = super(PlasmidEdit, self).get_form_kwargs()
        kwargs.update({'user': self.request.user})
        return kwargs

    def get_success_url(self, **kwargs):
        return reverse('plasmid', args=(self.object.id,)) + '?form_result_plasmid_edit_success=true'


class PlasmidCreate(CreateView):
    model = Plasmid
    form_class = PlasmidCreateForm
    template_name_suffix = '_create_form'

    @method_decorator(require_member_can_write_or_admin_current_project)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_form_kwargs(self):
        kwargs = super(PlasmidCreate, self).get_form_kwargs()
        kwargs.update({'user': self.request.user})
        return kwargs

    def form_valid(self, form):
        form.instance.project = get_current_project(self.request)
        return super().form_valid(form)

    def get_success_url(self, **kwargs):
        auto_create = ""
        if self.request.GET.get('b'):
            # comes from wizard --> auto assemble
            auto_create = "&ac"
        return reverse('plasmid', args=(self.object.id,)) + '?form_result_plasmid_create_success=true' + auto_create


class PlasmidCreateWizard(CreateView):
    model = Plasmid
    fields = '__all__'
    template_name_suffix = '_create_wizard'

    @method_decorator(require_member_can_write_or_admin_current_project)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["create_url"] = reverse('plasmid_create')
        context["next_url"] = reverse('plasmid_create_wizard_end')
        return context

    def form_valid(self, form):
        form.instance.project = get_current_project(self.request)
        return super().form_valid(form)

    def get_success_url(self, **kwargs):
        return reverse('plasmid', args=(self.object.id,))


@require_member_can_write_or_admin_current_project
def plasmid_create_wizard_end(request):
    context = {}
    if request.method == 'POST' and 'params' in request.POST:
        params = {}
        for param in request.POST.get('params').split('&'):
            name, values = param.split('=')
            if len(values.split('+')) > 1:
                values = values.split('+')
            params[name] = values
        if 'n' in params and 'i' in params and params['n'] and len(params['i']):
            backbone = None
            level = None
            sequencing_state = 0
            if 'b' in params:
                backbone = Plasmid.objects.get(id=params['b'])
                level = backbone.level
            else:
                # assume level from inserts
                first_insert = Plasmid.objects.get(id=params['i'][0])
                if first_insert.level is not None:
                    level = first_insert.level + 1
            if level:
                if level == 0:
                    sequencing_state = 1  # required
            description = ""
            if 'd' in params:
                description = params['d']
            plasmid_created = Plasmid.objects.create(
                name=params['n'],
                description=description,
                backbone=backbone,
                type=PlasmidType.objects.get(id=0),
                level=level,
                sequencing_state=sequencing_state,
                project=get_current_project(request),
            )
            for i in params['i']:
                plasmid_created.inserts.add(i)
            if backbone:
                for r in backbone.selectable_markers.all():
                    plasmid_created.selectable_markers.add(r.id)
            context['plasmid_create_result'] = plasmid_create_from_inserts(plasmid_created)
            return redirect('plasmid', plasmid_id=plasmid_created.id)
        else:
            context['wizard_error'] = 'Name & backbone & inserts are required fields.'
    else:
        context['wizard_error'] = "Plasmid can\'t be created. No parameters set."
    return render(request, 'inventory/plasmid.html', context)


@require_member_can_write_or_admin_current_project
def PlasmidCreateL0d(request):
    # check parameters
    if request.POST and request.POST.get('oh5-name') and request.POST.get('oh5-oh') and request.POST.get('oh3-name') and request.POST.get('oh3-oh') and request.POST.get('name') and request.POST.get('seq'):
        # take the backbone
        try:
            the_re = RestrictionEnzyme.objects.get(name="SapI")
            backbone = Plasmid.objects.get(name="pL0R-lacZ")

            if not the_re or not backbone:
                raise ObjectDoesNotExist

            plasmid_backbone_seq_result = seqio_get(backbone)
            oh_length = abs(the_re.rcut - the_re.fcut)

            if plasmid_backbone_seq_result[0] and the_re:
                backbone_record = plasmid_backbone_seq_result[1]
                hits = re_find_cut_positions(backbone_record.seq, the_re, True, True)
                if len(hits) == 2 and hits[1] > hits[0]:

                    final_record = backbone_record[0:hits[0] + oh_length - 1]

                    # go for the inserts
                    rec_oh_5 = SeqRecord(
                        Seq(request.POST.get('oh5-oh').upper()),
                        id=request.POST.get('oh5-name'),
                        annotations={"molecule_type": "DNA"}
                    )
                    rec_oh_5.features.append(SeqFeature(
                        FeatureLocation(0, len(request.POST.get('oh5-oh'))),
                        type="misc_feature",
                        strand=1,
                        qualifiers={
                            'ApEinfo_label': request.POST.get('oh5-name'),
                            'label': request.POST.get('oh5-name'),
                            'locus_tag': request.POST.get('oh5-name'),
                        }
                    ))
                    final_record = final_record + rec_oh_5


                    rec_seq = SeqRecord(
                        Seq(request.POST.get('seq').upper()),
                        id=request.POST.get('name'),
                        annotations={"molecule_type": "DNA"}
                    )
                    rec_seq.features.append(SeqFeature(
                        FeatureLocation(0, len(request.POST.get('seq'))),
                        type="misc_feature",
                        strand=1,
                        qualifiers={
                            'ApEinfo_label': request.POST.get('name'),
                            'label': request.POST.get('name'),
                            'locus_tag': request.POST.get('name'),
                        }
                    ))
                    final_record = final_record + rec_seq

                    if request.POST.get('oh3-prev_bases'):
                        prev_bases = request.POST.get('oh3-prev_bases')
                        rec_oh3_tc = SeqRecord(
                            Seq(prev_bases.upper()),
                            id=prev_bases.upper(),
                            annotations={"molecule_type": "DNA"}
                        )
                        rec_oh3_tc.features.append(SeqFeature(
                            FeatureLocation(0, len(prev_bases)),
                            type="misc_feature",
                            strand=1,
                            qualifiers={
                                'ApEinfo_label': prev_bases.upper(),
                                'label': prev_bases.upper(),
                                'locus_tag': prev_bases.upper(),
                            }
                        ))
                        final_record = final_record + rec_oh3_tc

                    if request.POST.get('oh3-stop'):
                        rec_oh3_stop = SeqRecord(
                            Seq(request.POST.get('oh3-stop').upper()),
                            id='STOP',
                            annotations={"molecule_type": "DNA"}
                        )
                        rec_oh3_stop.features.append(SeqFeature(
                            FeatureLocation(0, len(request.POST.get('oh3-stop'))),
                            type="misc_feature",
                            strand=1,
                            qualifiers={
                                'ApEinfo_label': 'STOP',
                                'label': 'STOP',
                                'locus_tag': 'STOP',
                            }
                        ))
                        final_record = final_record + rec_oh3_stop

                    rec_oh_3 = SeqRecord(
                        Seq(request.POST.get('oh3-oh').upper()),
                        id=request.POST.get('oh3-name'),
                        annotations={"molecule_type": "DNA"}
                    )
                    rec_oh_3.features.append(SeqFeature(
                        FeatureLocation(0, len(request.POST.get('oh3-oh'))),
                        type="misc_feature",
                        strand=1,
                        qualifiers={
                            'ApEinfo_label': request.POST.get('oh3-name'),
                            'label': request.POST.get('oh3-name'),
                            'locus_tag': request.POST.get('oh3-name'),
                        }
                    ))
                    final_record = final_record + rec_oh_3

                    final_record = final_record + backbone_record[hits[1] - 1:]
                    final_record.name = request.POST.get('name')
                    final_record.description = "Created with Weaver L0 Designer"
                    final_record.annotations = {"molecule_type": "DNA", "topology": "circular"}

                    plasmid_created = Plasmid.objects.create(
                        name=request.POST.get('name'),
                        description="Created with Weaver L0 Designer",
                        backbone=backbone,
                        type=PlasmidType.objects.get(id=0),
                        level=backbone.level,
                        sequencing_state=1,  # required
                        project=get_current_project(request),
                    )
                    for r in backbone.selectable_markers.all():
                        plasmid_created.selectable_markers.add(r.id)
                    plasmid_created.sequence.save(final_record.name.replace(" ", "_") + ".gb", ContentFile(final_record.format("gb")))

                    return HttpResponseRedirect(reverse('plasmid', args=(plasmid_created.id,)))
                else:
                    context = {
                        'error': 'Backbone plasmid does not contains appropriate restriction sites'
                    }

        except:
            context = {
                'error': 'Backbone plasmid or Retriction Enzyme not found'
            }
    else:
        context = {
            'error': 'Bad input parameters'
        }
    return render(request, 'inventory/plasmid_create_l0d.html', context)


@csrf_exempt
@require_member_can_read_project_of_plasmid
def plasmid_view_edit(request, plasmid_id):
    try:
        plasmid_to_detail = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    warnings = []

    if member_can_write_or_admin_plasmid(plasmid_to_detail, request.user):
        if request.method == 'POST' and 'saveOve' in request.POST:
            if 'gbContent' in request.POST:
                # saving from OVE
                plasmid_to_detail.sequence.save(plasmid_to_detail.name + ".gb", ContentFile(request.POST['gbContent']))
                result = 'File saved'
            else:
                result = 'Error: no gbContent'
            return HttpResponse(json.dumps({
                'result': result
            }), content_type="application/json")

        if request.method == 'POST' and 'create' in request.POST:
            if not plasmid_to_detail.sequence:
                plasmid_to_detail.sequence.save(plasmid_to_detail.name + ".gb", ContentFile(
                    '''LOCUS       ''' + plasmid_to_detail.name.replace(" ", "_") + '''                   0 bp    DNA     circular  31-AUG-2021
    ORIGIN      
    //'''))
            else:
                warnings.append('Can\'t create empty sequence on this plasmid, already has one')

    else:
        warnings.append('Current user can\'t make modifications')

    sequence = grab_seq(plasmid_to_detail)
    digest_picker_enzymes = []
    if sequence[0]:
        digest_picker_enzymes = enzymes_with_effective_cuts(
            str(sequence[1]),
            RestrictionEnzyme.objects.all().order_by('name'),
            is_circular=True,
        )

    context = {
        'plasmid': plasmid_to_detail,
        'sequence_file_contents': plasmid_sequence_file_contents(plasmid_to_detail),
        'warnings': warnings,
        'RESTRICTION_ENZYMES': digest_picker_enzymes,
        'user_can_edit_plasmid': member_can_write_or_admin_plasmid(plasmid_to_detail, request.user)
    }
    return render(request, 'inventory/plasmid_view_edit.html', context)


def plasmid_sequence_file_contents(plasmid):
    with open(plasmid.sequence.path, 'r') as file:
        sequence_file_contents = html.unescape(file.read())

    return re.sub(r'[\'"]', '', sequence_file_contents)


def plasmid_seqrecord(plasmid):
    if not plasmid.sequence:
        return None
    try:
        with open(plasmid.sequence.path, "r") as handle:
            records = list(SeqIO.parse(handle, "genbank"))
    except Exception:
        return None
    return records[0] if records else None


def parse_required_regions(raw_regions, sequence_length):
    zero_based_regions = [
        {
            'start': int(region['start']) - 1,
            'end': int(region['end']) - 1,
        }
        for region in raw_regions
        if str(region.get('start', '')).strip() and str(region.get('end', '')).strip()
    ]
    return normalize_regions(zero_based_regions, sequence_length)


def optional_int_query_param(request, name, default=None):
    raw_value = str(request.GET.get(name, '')).strip()
    return int(raw_value) if raw_value else default


def interval_contains(container_start, container_end, contained_start, contained_end, flank_bp=0):
    return container_start - flank_bp <= contained_start and container_end + flank_bp >= contained_end


def amplicon_contains_region(amplicon, region, sequence_length, flank_bp=0):
    amplicon_ranges = amplicon_segments(amplicon, sequence_length)
    region_ranges = (
        [(region.start, region.end)]
        if region.start <= region.end
        else [(region.start, sequence_length - 1), (0, region.end)]
    )
    for region_start, region_end in region_ranges:
        if not any(
                interval_contains(amplicon_start, amplicon_end, region_start, region_end, flank_bp)
                for amplicon_start, amplicon_end in amplicon_ranges):
            return False
    return True


def amplicon_matches_required_regions(amplicon, regions, sequence_length, flank_bp=0):
    return all(amplicon_contains_region(amplicon, region, sequence_length, flank_bp) for region in regions)


def amplicon_matches_primer_id(amplicon, primer_id):
    if not primer_id:
        return True
    notes = amplicon.get("notes") or {}
    return primer_id in (
        (notes.get("fwd_primer_id") or [""])[0],
        (notes.get("rev_primer_id") or [""])[0],
    )


def amplicon_matches_any_primer_id(amplicon, primer_ids):
    primer_ids = [str(primer_id).strip() for primer_id in primer_ids if str(primer_id).strip()]
    if not primer_ids:
        return True
    return any(amplicon_matches_primer_id(amplicon, primer_id) for primer_id in primer_ids)


@require_member_can_read_project_of_plasmid
def api_plasmid_primer_matches(request, plasmid_id):
    try:
        plasmid_to_match = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    sequence = grab_seq(plasmid_to_match)
    if not sequence[0]:
        return JsonResponse({
            'error': sequence[1],
            'primers': [],
            'count': 0
        }, status=400)

    primers = visible_primers_for_user(request.user)
    annotations = matching_primer_annotations(str(sequence[1]), primers)
    return JsonResponse({
        'primers': annotations,
        'count': len(annotations)
    })


@require_member_can_read_project_of_plasmid
def api_plasmid_amplicon_matches(request, plasmid_id):
    try:
        plasmid_to_match = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    sequence = grab_seq(plasmid_to_match)
    if not sequence[0]:
        return JsonResponse({
            'error': sequence[1],
            'amplicons': [],
            'count': 0
        }, status=400)

    try:
        min_product_size = optional_int_query_param(request, 'min_size', 100)
        max_product_size = optional_int_query_param(request, 'max_size')
        region_flank_bp = optional_int_query_param(
            request,
            'region_flank_bp',
            DEFAULT_AMPLICON_REGION_FLANK_BP,
        )
        max_tm_difference = float(request.GET.get('max_tm_diff', 5))
        raw_regions = json.loads(request.GET.get('regions', '[]'))
        required_regions = parse_required_regions(raw_regions, len(str(sequence[1])))
        raw_primer_ids = []
        for value in request.GET.getlist('primer_ids'):
            raw_primer_ids.extend(value.split(','))
        legacy_primer_id = str(request.GET.get('primer_id', '')).strip()
        if legacy_primer_id:
            raw_primer_ids.append(legacy_primer_id)
        primer_ids = tuple(dict.fromkeys(
            primer_id.strip()
            for primer_id in raw_primer_ids
            if primer_id.strip()
        ))
    except (TypeError, ValueError, KeyError, json.JSONDecodeError) as error:
        return JsonResponse({
            'error': 'Bad amplicon filter parameters: ' + str(error),
            'amplicons': [],
            'count': 0
        }, status=400)

    primers = visible_primers_for_user(request.user)
    candidate_annotations = matching_amplicon_annotations(
        str(sequence[1]),
        primers,
        min_product_size=min_product_size,
        max_product_size=max_product_size,
        max_tm_difference=max_tm_difference,
    )
    if required_regions:
        candidate_annotations = [
            amplicon for amplicon in candidate_annotations
            if amplicon_matches_required_regions(
                amplicon,
                required_regions,
                len(str(sequence[1])),
                region_flank_bp,
            )
        ]
    if primer_ids:
        candidate_annotations = [
            amplicon for amplicon in candidate_annotations
            if amplicon_matches_any_primer_id(amplicon, primer_ids)
        ]
    non_overlapping = request.GET.get('non_overlapping', 'true').lower() not in ('0', 'false', 'no')
    annotations = select_non_overlapping_amplicons(
        candidate_annotations,
        len(str(sequence[1]))
    ) if non_overlapping else candidate_annotations
    return JsonResponse({
        'amplicons': annotations,
        'candidates': candidate_annotations,
        'count': len(annotations),
        'candidate_count': len(candidate_annotations),
        'non_overlapping': non_overlapping,
        'filters': {
            'min_size': min_product_size,
            'max_size': max_product_size,
            'primer_ids': list(primer_ids),
            'region_count': len(required_regions),
            'region_flank_bp': region_flank_bp,
        },
    })


@require_member_can_read_project_of_plasmid
def api_plasmid_restriction_digests(request, plasmid_id):
    try:
        plasmid_to_digest = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    sequence = grab_seq(plasmid_to_digest)
    if not sequence[0]:
        return JsonResponse({
            'error': sequence[1],
            'results': [],
            'count': 0
        }, status=400)

    sequence_text = str(sequence[1])
    try:
        raw_required_enzymes = []
        for value in request.GET.getlist('required_enzymes'):
            raw_required_enzymes.extend(value.split(','))
        required_enzymes = tuple(dict.fromkeys(
            enzyme.strip()
            for enzyme in raw_required_enzymes
            if enzyme.strip()
        ))
        raw_regions = json.loads(request.GET.get('regions', '[]'))
        constraints = DigestConstraints(
            min_fragments=max(1, int(request.GET.get('min_fragments', DEFAULT_MIN_FRAGMENTS))),
            max_fragments=max(1, int(request.GET.get('max_fragments', DEFAULT_MAX_FRAGMENTS))),
            min_band_difference_bp=max(0, int(request.GET.get('min_band_difference_bp', DEFAULT_MIN_BAND_DIFFERENCE_BP))),
            min_fragment_size_bp=max(0, int(request.GET.get('min_fragment_size_bp', DEFAULT_MIN_FRAGMENT_SIZE_BP))),
            min_buffer_activity_percent=DEFAULT_MIN_BUFFER_ACTIVITY_PERCENT,
            max_enzymes=max(1, min(2, int(request.GET.get('max_enzymes', DEFAULT_MAX_ENZYMES)))),
            limit=max(1, min(50, int(request.GET.get('limit', DEFAULT_RESULT_LIMIT)))),
            required_regions=parse_required_regions(raw_regions, len(sequence_text)),
            required_enzymes=required_enzymes,
        )
    except (TypeError, ValueError, KeyError, json.JSONDecodeError) as error:
        return JsonResponse({
            'error': 'Bad restriction digest parameters: ' + str(error),
            'results': [],
            'count': 0
        }, status=400)

    enzymes = list(RestrictionEnzyme.objects.all())
    available_enzyme_names = {
        enzyme.name.lower(): enzyme.name
        for enzyme in enzymes
    }
    unknown_required_enzymes = [
        enzyme
        for enzyme in constraints.required_enzymes
        if enzyme.lower() not in available_enzyme_names
    ]
    if unknown_required_enzymes:
        return JsonResponse({
            'error': 'Bad restriction digest parameters: required enzyme not available: ' + ', '.join(unknown_required_enzymes),
            'results': [],
            'count': 0
        }, status=400)

    if constraints.max_fragments < constraints.min_fragments:
        return JsonResponse({
            'error': 'Bad restriction digest parameters: maximum fragments must be greater than or equal to minimum fragments',
            'results': [],
            'count': 0
        }, status=400)

    response = serialize_digest_response(
        sequence_text,
        enzymes,
        constraints=constraints,
        is_circular=True,
    )
    return JsonResponse(response)


@require_member_can_read_project_of_plasmid
def plasmid_download(request, plasmid_id):
    try:
        plasmid_to_download = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    if request.method == 'GET' and 'format' in request.GET:
        record_response = seqio_get(plasmid_to_download)
        if record_response[0]:
            response = HttpResponse(record_response[1].format(request.GET['format']), content_type="plain/text")
            response['Content-Disposition'] = 'inline; filename=' + plasmid_to_download.__str__() + '.' + request.GET[
                'format']
            return response

    # return original file if no format is specified
    file_path = os.path.join(settings.MEDIA_ROOT, plasmid_to_download.sequence.name)
    if os.path.exists(file_path):
        with open(file_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="plain/text")
            response['Content-Disposition'] = 'inline; filename=' + plasmid_to_download.__str__() + \
                                              os.path.splitext(os.path.basename(file_path))[1]
            return response
    raise Http404


@require_member_can_read_project_of_plasmid
def plasmid_download_clustal(request, plasmid_id):
    try:
        plasmid_to_download = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    file_path = os.path.join(settings.MEDIA_ROOT, plasmid_to_download.sequencing_clustal_file.name)
    if os.path.exists(file_path):
        with open(file_path, 'rb') as fh:
            response = HttpResponse(fh.read(), content_type="plain/text")
            response['Content-Disposition'] = 'inline; filename=' + plasmid_to_download.name + \
                                              os.path.splitext(os.path.basename(file_path))[1]
            return response
    raise Http404


@require_member_can_read_project_of_plasmid
def plasmid_label(request, plasmid_id):
    try:
        plasmid_to_label = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    context = {
        'plasmid': plasmid_to_label,
        'form': PlasmidLabel(),
        'user_can_edit_plasmid': member_can_write_or_admin_plasmid(plasmid_to_label, request.user)
    }
    return render(request, 'inventory/plasmid_label.html', context)


@require_member_can_read_project_of_plasmid
def plasmid_digest(request, plasmid_id):
    try:
        plasmid_to_digest = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    context = {
        'plasmid': plasmid_to_digest,
        'user_can_edit_plasmid': member_can_write_or_admin_plasmid(plasmid_to_digest, request.user)
    }

    sequence = grab_seq(plasmid_to_digest)

    if sequence[0]:
        context['digest_form'] = DigestForm()
        res = RestrictionEnzyme.objects.all()
        for the_re in res:
            the_re.hits, the_re.fragments = re_find_cut_fragments(sequence[1], the_re, True)

        context['res'] = res

        if request.method == 'POST':
            selected_res = []
            post_enzymes = json.loads(request.POST['enzymes'])
            for post_enzyme in post_enzymes:
                the_re = RestrictionEnzyme.objects.filter(name=post_enzyme).first()
                if the_re:
                    selected_res.append(the_re)
            context['selected_res'] = selected_res
            context['fragments'] = re_digestion_fragments(sequence[1], selected_res, True)

    else:
        context['error'] = sequence[1]

    return render(request, 'inventory/plasmid_digest.html', context)


@require_member_can_read_project_of_plasmid
def plasmid_pcr(request, plasmid_id):
    try:
        plasmid_to_pcr = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    context = {
        'plasmid': plasmid_to_pcr,
        'show_new_PCR': False,
        'user_can_edit_plasmid': member_can_write_or_admin_plasmid(plasmid_to_pcr, request.user)
    }

    sequence = grab_seq(plasmid_to_pcr)

    if sequence[0]:
        if request.method == 'GET' and 'start' in request.GET and 'end' in request.GET:
            try:
                selection_start = int(request.GET.get('start'))
                selection_end = int(request.GET.get('end'))
                margin = int(request.GET.get('margin', 300))
                max_tm_difference = float(request.GET.get('max_tm_diff', 5))
                if margin < 0:
                    margin = 0
            except ValueError:
                context['error'] = "Bad PCR design coordinates"
                return render(request, 'inventory/plasmid_pcr.html', context)

            primers = visible_primers_for_user(request.user)
            context['pcr_design'] = {
                'start': selection_start,
                'end': selection_end,
                'start_display': selection_start + 1,
                'end_display': selection_end + 1,
                'margin': margin,
                'max_tm_difference': max_tm_difference,
                'suggestions': suggest_pcr_primers(
                    str(sequence[1]),
                    primers,
                    selection_start,
                    selection_end,
                    margin=margin,
                    max_tm_difference=max_tm_difference,
                )
            }
            context['pcr_form'] = PCRForm(user=request.user)

        if request.method == 'POST':
            visible_primers = visible_primers_for_user(request.user)
            if request.POST['primer_f'] != "":
                primer_f = visible_primers.filter(id=request.POST['primer_f']).first()
                if primer_f is None:
                    context['error'] = "Forward primer is not available"
            else:
                if request.POST['primer_f_seq'] != "":
                    primer_f = Primer(
                        id='custom_f',
                        name='Custom F',
                        sequence_3=request.POST['primer_f_seq'],
                        sequence_5='',
                        fwd_or_rev='f',
                        intended_use='Custom sequence for PCR prediction'
                    )
                else:
                    context['error'] = "No forward primer set"
            if request.POST['primer_r'] != "":
                primer_r = visible_primers.filter(id=request.POST['primer_r']).first()
                if primer_r is None:
                    context['error'] = "Reverse primer is not available"
            else:
                if request.POST['primer_r_seq'] != "":
                    primer_r = Primer(
                        id='custom_r',
                        name='Custom R',
                        sequence_3=request.POST['primer_r_seq'],
                        sequence_5='',
                        fwd_or_rev='f',
                        intended_use='Custom sequence for PCR prediction'
                    )
                else:
                    context['error'] = "No forward primer set"
            if not 'error' in context:
                context['primer_f'] = primer_f
                context['primer_r'] = primer_r
                if primer_r.sequence_5:
                    context['primer_r_5_rc'] = str(Seq(primer_r.sequence_5).reverse_complement())
                if primer_r.sequence_3:
                    context['primer_r_3_rc'] = str(Seq(primer_r.sequence_3).reverse_complement())
                double_seq = str(sequence[1]) + str(sequence[1])
                pos_f = re.search(primer_f.sequence_3, double_seq, re.IGNORECASE)
                if pos_f:
                    start = pos_f.end()
                    seq_from_f = double_seq[start:]
                    pos_r = re.search(context['primer_r_3_rc'], seq_from_f, re.IGNORECASE)
                    if pos_r:
                        end = pos_r.start()
                        context['amplicon'] = seq_from_f[:end].lower()
                        context['size'] = len(
                            primer_f.sequence_5 + primer_f.sequence_3 + context['amplicon'] + primer_r.sequence_3 +
                            primer_r.sequence_5)
                    else:
                        context['error'] = "REV primer does not hit template"
                else:
                    context['error'] = "FWD primer does not hit template"
            context['show_new_PCR'] = True
        elif 'pcr_form' not in context:
            context['pcr_form'] = PCRForm(user=request.user)

    else:
        context['error'] = sequence[1]

    return render(request, 'inventory/plasmid_pcr.html', context)


def plasmid_save_clustal(plasmid, records):
    try:
        file_name = os.path.join("uploads", "sequencing_clustal", plasmid.name + ".clustal")
        output_handle = StringIO()
        align = Bio.Align.MultipleSeqAlignment(records)
        AlignIO.write(align, output_handle, "clustal")
        plasmid.sequencing_clustal_file.save(file_name, ContentFile(output_handle.getvalue()), save=True)
        return True, 'Save clustal file done'
    except Exception as e:
        return False, e.__str__()


def get_optimal_alignment(ref_seq, query_seq, is_reversed=False):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")

    if is_reversed:
        query_seq = reverse_complement(query_seq)

    try:
        alignments = aligner.align(ref_seq.upper(), query_seq.upper())
        return True, next(alignments)
    except MemoryError:
        return False, 'Too many alignments.'
    except Exception as e:
        return False, e.__str__()


def is_fasta_alignment_file(filename):
    return filename.lower().endswith((".fa", ".fas", ".fasta"))


def fasta_records_from_text(text):
    sequence_text = str(text or '').strip()
    if not sequence_text:
        return []
    fasta_text = sequence_text if ">" in sequence_text else f">Query\n{sequence_text}\n"
    fasta_text = "\n".join(
        line for line in fasta_text.splitlines()
        if line.strip() and not line.lstrip().startswith((";", "#", "!"))
    )
    fasta_io = StringIO(fasta_text)
    try:
        return list(SeqIO.parse(fasta_io, "fasta"))
    except ValueError:
        return []
    finally:
        fasta_io.close()


def render_uploaded_fasta_alignment(request, plasmid_to_align, upload_files, form, context):
    records = []
    for uploaded_file in upload_files:
        try:
            fasta_text = uploaded_file.read().decode("utf-8-sig")
        except UnicodeDecodeError:
            context['error'] = "{} is not valid UTF-8 text".format(uploaded_file.name)
            return render(request, 'inventory/plasmid_align_sanger.html', context)
        records.extend(fasta_records_from_text(fasta_text))
    if not records:
        context['error'] = "Input sequence not in FASTA format"
        return render(request, 'inventory/plasmid_align_sanger.html', context)
    plasmid_seq = grab_seq(plasmid_to_align)[1]
    fasta_result = fasta_alignment_result(str(plasmid_seq), records)
    if not any(read.get("is_usable") for read in fasta_result["reads"]):
        context['error'] = "No FASTA sequences aligned"
        return render(request, 'inventory/plasmid_align_sanger.html', context)

    fasta_result["reference_record"] = plasmid_seqrecord(plasmid_to_align)
    context['sanger_result'] = fasta_result
    context['sanger_browser_data'] = json.dumps(sanger_browser_data(str(plasmid_seq), fasta_result, plasmid_to_align.name))
    context['align_data'] = alignment_tracks_for_ove(plasmid_to_align.name, str(plasmid_seq), fasta_result["reads"])
    context['plasmid_sequence_file_contents'] = plasmid_sequence_file_contents(plasmid_to_align)
    context['alignment_source_type'] = "FASTA"
    context['fasta_view_mode'] = form.cleaned_data.get("alignment_view_mode", "combined")
    context['fasta_clustal_content'] = clustal_content(plasmid_to_align.name, str(plasmid_seq), fasta_result["reads"])
    context['fasta_clustal_filename'] = get_valid_filename("{}-fasta-alignment.clustal".format(plasmid_to_align.name))
    context['show_results'] = True
    return render(request, 'inventory/plasmid_align_sanger.html', context)


def fasta_alignment_result(reference_sequence, records):
    params = SangerProcessingParameters()
    if isinstance(records, SeqRecord):
        records = [records]
    reads = []
    for index, record in enumerate(records or [], start=1):
        read_sequence = str(record.seq).upper()
        qualities = [40] * len(read_sequence)
        alignment = align_read(
            reference_sequence,
            read_sequence,
            qualities,
            0,
            params,
            len(read_sequence),
        )
        record_name = record.id or record.name or "FASTA {}".format(index)
        reads.append({
            "name": record_name,
            "files": [],
            "formats": ["FASTA"],
            "parsed_sources": [],
            "selected_source": "fasta",
            "raw_sequence": read_sequence,
            "trimmed_sequence": read_sequence,
            "trim_start": 0,
            "trim_end": len(read_sequence),
            "quality_metrics": {
                "raw_length": len(read_sequence),
                "trimmed_length": len(read_sequence),
                "quality_available": False,
                "mean_quality": None,
                "alignment_blocks": [],
                "low_confidence_regions": [],
                "intermediate_confidence_regions": [],
            },
            "alignment": alignment,
            "warnings": [],
            "errors": [],
            "is_usable": bool(alignment),
            "unusable_reason": "" if alignment else "FASTA sequence did not align",
            "chromatogram": {},
        })
    combined = combined_metrics(reference_sequence, reads, params)
    classification = classify_run(combined, reads, params)
    if classification.get("reasons") == ["The high-quality Sanger-aligned region is consistent with the expected plasmid sequence"]:
        classification["reasons"] = ["The aligned FASTA sequence(s) are consistent with the expected plasmid sequence"]
    return {
        "parameters": params.as_dict(),
        "uploaded_files": [],
        "reads": reads,
        "combined": combined,
        "classification": classification,
    }


SANGER_FEATURE_TYPE_COLORS = {
    "cds": "#2fb344",
    "promoter": "#f0b429",
    "terminator": "#d94841",
    "restriction_site": "#8f63d9",
    "cut_site": "#8f63d9",
    "recombination_site": "#e66fb2",
    "primer_bind": "#2f9ed8",
    "rep_origin": "#20a39e",
    "origin_of_replication": "#20a39e",
}


def sanger_feature_role(feature_type, label):
    normalized_type = (feature_type or "").strip().lower()
    normalized_label = (label or "").strip().lower()
    enzyme_label = re.sub(r"[-_ ]?hf$", "", (label or "").strip(), flags=re.IGNORECASE)
    if normalized_type == "cds":
        return "cds"
    if "promoter" in normalized_type or "promoter" in normalized_label or normalized_label.startswith(("p_", "p-")):
        return "promoter"
    if "terminator" in normalized_type or "terminator" in normalized_label or normalized_label.startswith(("t_", "t-")):
        return "terminator"
    if (
        normalized_type in ("restriction_site", "cut_site")
        or "restriction" in normalized_type
        or "cut site" in normalized_label
        or enzyme_label in rest_dict
    ):
        return "restriction_site"
    if (
        normalized_type in ("recombination_site", "misc_recomb")
        or "recombination" in normalized_type
        or "recombination" in normalized_label
        or re.search(r"\bhr\d*[-_ ]?chr\d+\b", normalized_label)
        or normalized_label.startswith(("attb", "attp", "attl", "attr"))
        or normalized_label in ("loxp", "frt")
    ):
        return "recombination_site"
    if normalized_type in ("primer_bind", "primer") or "primer" in normalized_label:
        return "primer_bind"
    if normalized_type in ("rep_origin", "origin_of_replication", "ori") or normalized_label in ("ori", "origin") or "origin of replication" in normalized_label:
        return "rep_origin"
    return normalized_type


def sanger_feature_color(feature_type, label, qualifiers):
    role = sanger_feature_role(feature_type, label)
    if role in SANGER_FEATURE_TYPE_COLORS:
        return SANGER_FEATURE_TYPE_COLORS[role]
    explicit_color = (
        qualifiers.get("ApEinfo_fwdcolor", [""])[0]
        or qualifiers.get("ApEinfo_revcolor", [""])[0]
        or qualifiers.get("color", [""])[0]
    )
    return explicit_color or "#7fb3ff"


def sanger_browser_data(reference_sequence, service_result, reference_name="Reference"):
    def feature_rows():
        features = []
        record = service_result.get("reference_record")
        if not record:
            return features
        for index, feature in enumerate(getattr(record, "features", []) or []):
            if feature.type == "source":
                continue
            parts = list(feature.location.parts) if isinstance(feature.location, CompoundLocation) else [feature.location]
            label = (
                feature.qualifiers.get("label", [""])[0]
                or feature.qualifiers.get("gene", [""])[0]
                or feature.qualifiers.get("product", [""])[0]
                or feature.type
            )
            role = sanger_feature_role(feature.type, label)
            color = sanger_feature_color(feature.type, label, feature.qualifiers)
            strand = feature.location.strand or 0
            for part in parts:
                start = int(part.start) % len(reference_sequence) if reference_sequence else 0
                end = (int(part.end) - 1) % len(reference_sequence) if reference_sequence else 0
                features.append({
                    "id": "{}-{}-{}".format(index, start, end),
                    "name": label,
                    "type": feature.type,
                    "role": role,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "color": color,
                    "crosses_origin": int(part.end) > len(reference_sequence) or start > end,
                })
        return features

    def trace_reference_points(chromatogram, projection_base_indices):
        base_positions = chromatogram.get("basePos") or []
        qualities = chromatogram.get("qualNums") or []
        points = []
        for coord, base_index in enumerate(projection_base_indices or []):
            if base_index is None:
                continue
            if 0 <= int(base_index) < len(base_positions):
                quality = qualities[int(base_index)] if int(base_index) < len(qualities) else None
                points.append({
                    "coordinate": coord,
                    "baseIndex": int(base_index),
                    "tracePos": base_positions[int(base_index)],
                    "quality": quality,
                })
        return points

    def trace_max_signal(chromatogram):
        values = []
        for key in ("aTrace", "cTrace", "gTrace", "tTrace"):
            values.extend(chromatogram.get(key) or [])
        return max(values) if values else 0

    def fallback_projection_base_indices(read, alignment):
        projection = alignment.get("reference_projection", "")
        if not projection or not reference_sequence:
            return []
        existing = alignment.get("reference_projection_base_indices")
        if existing and len(existing) == len(reference_sequence):
            return existing

        quality = read.get("quality_metrics") or {}
        trim_start = int(quality.get("trim_start") or 0)
        trim_end = int(quality.get("trim_end") or (trim_start + len(read.get("trimmed_sequence", ""))))
        chromatogram = read.get("chromatogram") or {}
        qualities = chromatogram.get("qualNums") or []
        if read.get("trimmed_sequence") and qualities:
            recalculated = align_read(
                reference_sequence,
                read.get("trimmed_sequence", ""),
                qualities,
                trim_start,
                SangerProcessingParameters(),
                trim_end,
            )
            if recalculated and recalculated.get("reference_projection_base_indices"):
                alignment.update({
                    "reference_projection": recalculated.get("reference_projection", projection),
                    "reference_projection_base_indices": recalculated.get("reference_projection_base_indices", []),
                    "query_start": recalculated.get("query_start"),
                    "query_end": recalculated.get("query_end"),
                    "best_orientation": recalculated.get("best_orientation"),
                })
                return recalculated["reference_projection_base_indices"]
        start = int(alignment.get("start") or 0)
        end = int(alignment.get("end") or 0)
        orientation = alignment.get("best_orientation") or alignment.get("orientation") or "forward"
        indices = [None] * len(reference_sequence)
        coords = []
        coord = start
        while True:
            coords.append(coord)
            if coord == end:
                break
            coord = (coord + 1) % len(reference_sequence)
            if len(coords) > len(reference_sequence):
                break
        query_index = 0
        for coord in coords:
            if coord >= len(projection) or projection[coord] == "-":
                continue
            indices[coord] = trim_end - 1 - query_index if orientation == "reverse" else trim_start + query_index
            query_index += 1
        return indices

    reads = []
    for read in service_result.get("reads", []):
        alignment = read.get("alignment") or {}
        display_alignment = read.get("display_alignment") or alignment
        if not read.get("is_usable") or not display_alignment:
            continue
        projection_base_indices = fallback_projection_base_indices(read, display_alignment)
        reads.append({
            "name": read["name"],
            "orientation": display_alignment.get("orientation", "unmapped"),
            "start": display_alignment.get("start", 0),
            "end": display_alignment.get("end", 0),
            "segments": display_alignment.get("segments", []),
            "crosses_origin": display_alignment.get("crosses_origin", False),
            "identity": alignment.get("identity", display_alignment.get("identity", 0)),
            "projection": display_alignment.get("reference_projection", ""),
            "projectionBaseIndices": projection_base_indices,
            "chromatogram": read.get("chromatogram") or {},
            "confidenceRegions": (read.get("quality_metrics") or {}).get("low_confidence_regions", []),
            "intermediateConfidenceRegions": (read.get("quality_metrics") or {}).get("intermediate_confidence_regions", []),
            "acceptedBlocks": (read.get("quality_metrics") or {}).get("alignment_blocks", []),
            "traceReferencePoints": trace_reference_points(read.get("chromatogram") or {}, projection_base_indices),
            "traceMaxSignal": trace_max_signal(read.get("chromatogram") or {}),
            "variants": [
                {
                    "coordinate": variant.get("coordinate", 0),
                    "type": variant.get("type", ""),
                    "expected": variant.get("expected", ""),
                    "observed": variant.get("observed", ""),
                    "quality": variant.get("quality"),
                    "low_quality": variant.get("low_quality", False),
                    "base_index": variant.get("base_index"),
                }
                for variant in display_alignment.get("variants", [])
            ],
        })
    return {
        "referenceName": reference_name,
        "referenceSequence": reference_sequence,
        "referenceLength": len(reference_sequence),
        "displayOrigin": reads[0]["start"] if reads else 0,
        "features": feature_rows(),
        "depth": service_result.get("combined", {}).get("depth", []),
        "reads": reads,
        "uncoveredRegions": service_result.get("combined", {}).get("uncovered_regions", []),
    }


def chromatogram_for_saved_read(read):
    stored = read.parsing_result.get("chromatogram", {})
    if read.selected_source != "ab1":
        return stored
    ab1_file = next((file_obj for file_obj in read.files.all() if file_obj.format == "ab1" and file_obj.file), None)
    if not ab1_file:
        return stored
    try:
        with ab1_file.file.open("rb") as handle:
            parsed = parse_ab1(handle.read())
    except Exception:
        return stored
    return parsed.chromatogram or stored


def recalculated_saved_sanger_read(read, reference_sequence, params=None):
    params = params or SangerProcessingParameters()
    chromatogram = chromatogram_for_saved_read(read)
    raw_sequence = read.raw_sequence or "".join(chromatogram.get("baseCalls") or [])
    qualities = chromatogram.get("qualNums") or []
    if not raw_sequence or not qualities:
        return None
    try:
        trimmed, trim_start, trim_end, quality_metrics = trim_by_quality(raw_sequence, qualities, params, chromatogram)
        usable, unusable_reason = read_is_usable(trimmed, quality_metrics, read.parsing_result.get("errors", []), params)
        alignment = align_read(
            reference_sequence,
            trimmed,
            qualities,
            trim_start,
            params,
            trim_end,
            trusted_blocks=quality_metrics.get("alignment_blocks"),
            forced_orientation=read.detected_orientation if read.detected_orientation in ("forward", "reverse") else None,
        ) if usable else None
    except Exception:
        return None
    display_alignment = None
    if alignment:
        display_start, display_end = display_trim_range(len(raw_sequence), quality_metrics)
        display_sequence = raw_sequence[display_start:display_end]
        if display_sequence:
            display_alignment = align_read(
                reference_sequence,
                display_sequence,
                qualities,
                display_start,
                params,
                display_end,
                forced_orientation=alignment.get("best_orientation") or alignment.get("orientation"),
            )
    warnings = list(read.warnings or [])
    if unusable_reason and unusable_reason not in warnings:
        warnings.append(unusable_reason)
    return {
        "id": str(read.id),
        "name": read.name,
        "formats": read.parsing_result.get("formats", []),
        "selected_source": read.selected_source,
        "raw_sequence": raw_sequence,
        "trimmed_sequence": trimmed,
        "quality_metrics": quality_metrics,
        "alignment": alignment or {},
        "display_alignment": display_alignment,
        "warnings": warnings,
        "errors": read.parsing_result.get("errors", []),
        "is_usable": usable and bool(alignment),
        "chromatogram": chromatogram,
    }


def sanger_result_from_run(run):
    reference_sequence = str(grab_seq(run.plasmid)[1])
    params = SangerProcessingParameters()

    variants = []
    variants_by_read = {}
    for variant in run.variants.select_related("read"):
        row = {
            "read": variant.read.name if variant.read else "",
            "coordinate": variant.coordinate,
            "type": variant.variant_type,
            "expected": variant.expected_base,
            "observed": variant.observed_base,
            "quality": variant.quality,
            "low_quality": "low_quality" in variant.flags,
            "base_index": (variant.evidence or {}).get("base_index"),
        }
        variants.append(row)
        if variant.read_id:
            variants_by_read.setdefault(variant.read_id, []).append(row)

    reads = []
    for read in run.reads.prefetch_related("files").all():
        chromatogram = chromatogram_for_saved_read(read)
        recalculated_read = recalculated_saved_sanger_read(read, reference_sequence, params)
        if recalculated_read:
            reads.append(recalculated_read)
            continue
        alignment = read.alignment_metrics.copy() if read.alignment_metrics else {}
        saved_variants = variants_by_read.get(read.id, [])
        alignment_variant_base_indices = {}
        for variant in alignment.get("variants", []):
            key = (
                variant.get("coordinate"),
                variant.get("type"),
                variant.get("observed", ""),
                variant.get("expected", ""),
            )
            if variant.get("base_index") is not None:
                alignment_variant_base_indices[key] = variant.get("base_index")
        for variant in saved_variants:
            if variant.get("base_index") is not None:
                continue
            key = (
                variant.get("coordinate"),
                variant.get("type"),
                variant.get("observed", ""),
                variant.get("expected", ""),
            )
            if key in alignment_variant_base_indices:
                variant["base_index"] = alignment_variant_base_indices[key]
        alignment["variants"] = saved_variants
        reads.append({
            "id": str(read.id),
            "name": read.name,
            "formats": read.parsing_result.get("formats", []),
            "selected_source": read.selected_source,
            "raw_sequence": read.raw_sequence,
            "trimmed_sequence": read.trimmed_sequence,
            "quality_metrics": read.quality_metrics,
            "alignment": alignment,
            "warnings": read.warnings,
            "errors": read.parsing_result.get("errors", []),
            "is_usable": read.is_usable,
            "chromatogram": chromatogram,
        })

    if any(read.get("alignment") for read in reads):
        combined = combined_metrics(reference_sequence, reads, params)
        variants = combined.get("variants", variants)
    else:
        combined = run.combined_metrics.copy() if run.combined_metrics else {}
        combined["variants"] = variants
    labels = {
        "PASS": "Verifica",
        "REVIEW": "Requiere revisión",
        "FAIL": "No verifica",
        "NO_DATA": "Sin datos utilizables",
    }
    return {
        "parameters": run.parameters,
        "reads": reads,
        "combined": combined,
        "classification": {
            "state": run.automated_state,
            "label": labels.get(run.automated_state, run.automated_state),
            "reasons": run.automated_reasons,
        },
    }


def persist_sanger_verification(plasmid, user, service_result, label="", notes="", save_clustal=False):
    run = SangerVerificationRun.objects.create(
        plasmid=plasmid,
        created_by=user,
        label=label or "",
        notes=notes or "",
        parameters=service_result["parameters"],
        automated_state=service_result["classification"]["state"],
        automated_reasons=service_result["classification"]["reasons"],
        combined_metrics={key: value for key, value in service_result["combined"].items() if key != "variants"},
    )
    saved_reads = {}
    for read_data in service_result["reads"]:
        alignment = read_data.get("alignment") or {}
        read = SangerRead.objects.create(
            run=run,
            name=read_data["name"],
            detected_orientation=alignment.get("orientation", "unmapped"),
            raw_sequence=read_data.get("raw_sequence", ""),
            trimmed_sequence=read_data.get("trimmed_sequence", ""),
            selected_source=read_data.get("selected_source", ""),
            parsing_result={
                "formats": read_data.get("formats", []),
                "errors": read_data.get("errors", []),
                "unusable_reason": read_data.get("unusable_reason", ""),
                "chromatogram": read_data.get("chromatogram", {}),
            },
            quality_metrics=read_data.get("quality_metrics", {}),
            alignment_metrics={key: value for key, value in alignment.items() if key not in ("covered_positions", "variants", "ref_alignment", "read_alignment", "oriented_sequence")},
            warnings=read_data.get("warnings", []),
            is_usable=read_data.get("is_usable", False),
        )
        saved_reads[read_data["name"]] = read
        for uploaded in read_data.get("files", []):
            read_file = SangerReadFile.objects.create(
                read=read,
                format=uploaded.format,
                original_name=uploaded.original_name,
                sha256=uploaded.sha256,
                size=uploaded.size,
                metadata={},
            )
            read_file.file.save(uploaded.original_name, ContentFile(uploaded.data), save=True)

    for variant in service_result["combined"].get("variants", []):
        read = saved_reads.get(variant.get("read"))
        SangerVariant.objects.create(
            run=run,
            read=read,
            coordinate=variant.get("coordinate", 0),
            variant_type=variant.get("type", ""),
            expected_base=variant.get("expected", ""),
            observed_base=variant.get("observed", ""),
            quality=variant.get("quality"),
            evidence={"read": variant.get("read", ""), "base_index": variant.get("base_index")},
            flags=["low_quality"] if variant.get("low_quality") else [],
        )

    if save_clustal:
        file_text = clustal_content(plasmid.name, grab_seq(plasmid)[1], service_result["reads"])
        file_name = "{}-{}-sanger.clustal".format(plasmid.name, run.id)
        run.clustal_file.save(file_name, ContentFile(file_text), save=True)
        plasmid.sequencing_clustal_file = run.clustal_file.name
        plasmid.save()
    return run


def delete_sanger_run_files(run):
    plasmid = run.plasmid
    clustal_name = run.clustal_file.name if run.clustal_file else ""
    if clustal_name:
        run.clustal_file.storage.delete(clustal_name)
        if plasmid.sequencing_clustal_file and plasmid.sequencing_clustal_file.name == clustal_name:
            plasmid.sequencing_clustal_file = None
            plasmid.save()
    for read_file in SangerReadFile.objects.filter(read__run=run):
        if read_file.file:
            read_file.file.storage.delete(read_file.file.name)


@require_member_can_read_project_of_plasmid
def plasmid_align_fasta(request, plasmid_id):
    try:
        plasmid_to_align = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    context = {
        'plasmid': plasmid_to_align,
        'user_can_edit_plasmid': member_can_write_or_admin_plasmid(plasmid_to_align, request.user),
        'alignment_source_type': "FASTA",
    }

    if request.method == 'POST':
        form = FastaAlignForm(request.POST, request.FILES)
        if form.is_valid():
            context['fasta_view_mode'] = form.cleaned_data.get("alignment_view_mode", "combined")
            records = []
            saw_input = False
            if request.POST.get('fasta_sequence'):
                saw_input = True
                records.extend(fasta_records_from_text(request.POST.get('fasta_sequence')))

            if request.FILES.getlist("fasta_file"):
                saw_input = True
                for fasta_file in request.FILES.getlist("fasta_file"):
                    try:
                        fasta_text = b"".join(fasta_file.chunks()).decode("utf-8-sig")
                    except UnicodeDecodeError:
                        context['error'] = '{} is not valid UTF-8 text'.format(fasta_file.name)
                        break
                    records.extend(fasta_records_from_text(fasta_text))

            if not records and not context.get('error'):
                context['error'] = 'Input sequence not in FASTA format' if saw_input else 'No input sequence'

            if records:
                plasmid_seq = grab_seq(plasmid_to_align)[1]
                fasta_result = fasta_alignment_result(str(plasmid_seq), records)

                if any(read.get("is_usable") for read in fasta_result["reads"]):
                    fasta_result["reference_record"] = plasmid_seqrecord(plasmid_to_align)
                    context['sanger_result'] = fasta_result
                    context['sanger_browser_data'] = json.dumps(sanger_browser_data(str(plasmid_seq), fasta_result, plasmid_to_align.name))
                    context['align_data'] = alignment_tracks_for_ove(plasmid_to_align.name, str(plasmid_seq), fasta_result["reads"])
                    context['alignment_source_type'] = "FASTA"
                    context['fasta_clustal_content'] = clustal_content(plasmid_to_align.name, str(plasmid_seq), fasta_result["reads"])
                    context['fasta_clustal_filename'] = get_valid_filename("{}-fasta-alignment.clustal".format(plasmid_to_align.name))
                    context['show_results'] = True

                    context['plasmid_sequence_file_contents'] = plasmid_sequence_file_contents(plasmid_to_align)
                else:
                    context['error'] = "No FASTA sequences aligned"
            else:
                if not context.get('error'):
                    context['error'] = 'Error while parsing input sequence'
    else:
        context['upload_form'] = FastaAlignForm()
        context['show_upload_form'] = True
    return render(request, 'inventory/plasmid_align_sanger.html', context)


@require_member_can_read_project_of_plasmid
def plasmid_align_sanger(request, plasmid_id):
    try:
        plasmid_to_align = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    context = {
        'plasmid': plasmid_to_align,
        'user_can_edit_plasmid': member_can_write_or_admin_plasmid(plasmid_to_align, request.user),
        'recent_runs': plasmid_to_align.sanger_verification_runs.all()[:10],
    }

    if request.method == 'POST':
        form = SangerAlignForm(request.POST, request.FILES)
        if form.is_valid():
            upload_files = request.FILES.getlist("sanger_files")
            legacy_ab1 = request.FILES.get("ab1")
            if legacy_ab1 and legacy_ab1 not in upload_files:
                upload_files.append(legacy_ab1)
            fasta_files = [uploaded for uploaded in upload_files if is_fasta_alignment_file(uploaded.name)]
            if fasta_files:
                if len(fasta_files) != len(upload_files):
                    context['error'] = "Do not mix FASTA files with Sanger trace files in the same alignment batch."
                    context['upload_form'] = form
                    context['show_upload_form'] = True
                    return render(request, 'inventory/plasmid_align_sanger.html', context)
                return render_uploaded_fasta_alignment(request, plasmid_to_align, upload_files, form, context)
            try:
                plasmid_seq = str(grab_seq(plasmid_to_align)[1])
                sanger_result = process_sanger_files(upload_files, plasmid_seq)
                sanger_result["reference_record"] = plasmid_seqrecord(plasmid_to_align)
                run = persist_sanger_verification(
                    plasmid_to_align,
                    request.user,
                    sanger_result,
                    label=form.cleaned_data.get("label", ""),
                    notes=form.cleaned_data.get("notes", ""),
                    save_clustal=form.cleaned_data.get("save_clustal_file", False),
                )
                context['run'] = run
                context['sanger_result'] = sanger_result
                context['align_data'] = alignment_tracks_for_ove(plasmid_to_align.name, plasmid_seq, sanger_result["reads"])
                context['sanger_browser_data'] = json.dumps(sanger_browser_data(plasmid_seq, sanger_result, plasmid_to_align.name))
                context['plasmid_sequence_file_contents'] = plasmid_sequence_file_contents(plasmid_to_align)
                context['show_results'] = True
                context['recent_runs'] = plasmid_to_align.sanger_verification_runs.exclude(id=run.id)[:10]
                if form.cleaned_data.get("save_clustal_file"):
                    context['save_clustal_done'] = "Saved Sanger Clustal file for this verification run"
            except Exception as exc:
                context['error'] = "Sanger verification failed: {}".format(exc)
        else:
            context['upload_form'] = form
            context['show_upload_form'] = True
    else:
        context['upload_form'] = SangerAlignForm()
        context['show_upload_form'] = True
    return render(request, 'inventory/plasmid_align_sanger.html', context)


def redirect_to_sanger_verification(request, plasmid_to_align):
    if not on_project_member_can_any(plasmid_to_align.project, request.user):
        return render(request, 'common/no_permission_to_edit.html')

    runs = plasmid_to_align.sanger_verification_runs.all()
    run = runs.filter(manual_decision="VERIFIED").first() or runs.first()
    if run:
        return redirect("sanger_run_detail", plasmid_id=plasmid_to_align.id, run_id=run.id)
    return redirect("plasmid_align_sanger", plasmid_id=plasmid_to_align.id)


def plasmid_seq_verification_entry(request, weaver_id):
    try:
        plasmid_to_align = Plasmid.objects.get(idx=weaver_id)
    except ObjectDoesNotExist:
        raise Http404
    return redirect_to_sanger_verification(request, plasmid_to_align)


def plasmid_seq_verification_entry_by_uuid(request, plasmid_id):
    try:
        plasmid_to_align = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404
    return redirect_to_sanger_verification(request, plasmid_to_align)


@require_member_can_read_project_of_plasmid
def sanger_run_detail(request, plasmid_id, run_id):
    try:
        plasmid_to_align = Plasmid.objects.get(id=plasmid_id)
        run = SangerVerificationRun.objects.get(id=run_id, plasmid=plasmid_to_align)
    except ObjectDoesNotExist:
        raise Http404

    plasmid_seq = str(grab_seq(plasmid_to_align)[1])
    sanger_result = sanger_result_from_run(run)
    sanger_result["reference_record"] = plasmid_seqrecord(plasmid_to_align)
    context = {
        'plasmid': plasmid_to_align,
        'user_can_edit_plasmid': member_can_write_or_admin_plasmid(plasmid_to_align, request.user),
        'recent_runs': plasmid_to_align.sanger_verification_runs.exclude(id=run.id)[:10],
        'run': run,
        'sanger_result': sanger_result,
        'sanger_browser_data': json.dumps(sanger_browser_data(plasmid_seq, sanger_result, plasmid_to_align.name)),
        'align_data': alignment_tracks_for_ove(plasmid_to_align.name, plasmid_seq, sanger_result["reads"]),
        'plasmid_sequence_file_contents': plasmid_sequence_file_contents(plasmid_to_align),
        'show_results': True,
        'is_saved_run_view': True,
    }
    return render(request, 'inventory/plasmid_align_sanger.html', context)


@require_member_can_read_project_of_plasmid
def sanger_read_chromatogram(request, plasmid_id, run_id, read_id):
    try:
        plasmid_to_align = Plasmid.objects.get(id=plasmid_id)
        run = SangerVerificationRun.objects.get(id=run_id, plasmid=plasmid_to_align)
        read = SangerRead.objects.get(id=read_id, run=run)
    except ObjectDoesNotExist:
        raise Http404

    plasmid_seq = str(grab_seq(plasmid_to_align)[1])
    read_data = recalculated_saved_sanger_read(read, plasmid_seq) or {
        "id": str(read.id),
        "name": read.name,
        "selected_source": read.selected_source,
        "quality_metrics": read.quality_metrics,
        "alignment": read.alignment_metrics,
        "warnings": read.warnings,
        "errors": read.parsing_result.get("errors", []),
        "is_usable": read.is_usable,
        "chromatogram": chromatogram_for_saved_read(read),
    }
    chromatogram = read_data.get("chromatogram") or {}
    if not chromatogram.get("aTrace"):
        raise Http404
    source_file = read.files.filter(format=read.selected_source).first() or read.files.first()
    context = {
        "plasmid": plasmid_to_align,
        "run": run,
        "read": read,
        "read_data": read_data,
        "source_file_name": source_file.original_name if source_file else read.name,
        "chromatogram_data": json.dumps({
            "readName": read.name,
            "sourceFileName": source_file.original_name if source_file else read.name,
            "orientation": (read_data.get("alignment") or {}).get("orientation", read.detected_orientation),
            "chromatogram": chromatogram,
            "qualityMetrics": read_data.get("quality_metrics") or {},
        }),
    }
    return render(request, "inventory/sanger_chromatogram.html", context)


@require_member_can_read_project_of_plasmid
def sanger_run_download(request, plasmid_id, run_id, kind):
    try:
        run = SangerVerificationRun.objects.get(id=run_id, plasmid_id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404
    if kind == "variants":
        rows = []
        for variant in run.variants.select_related("read"):
            rows.append({
                "read": variant.read.name if variant.read else "",
                "coordinate": variant.coordinate,
                "type": variant.variant_type,
                "expected": variant.expected_base,
                "observed": variant.observed_base,
                "quality": variant.quality,
                "low_quality": "low_quality" in variant.flags,
            })
        response = HttpResponse(variants_csv(rows), content_type="text/csv")
        response['Content-Disposition'] = 'attachment; filename="{}-sanger-variants.csv"'.format(run.id)
        return response
    if kind == "reads":
        reads = []
        for read in run.reads.all():
            reads.append({
                "name": read.name,
                "formats": read.parsing_result.get("formats", []),
                "selected_source": read.selected_source,
                "is_usable": read.is_usable,
                "alignment": read.alignment_metrics,
                "quality_metrics": read.quality_metrics,
                "warnings": read.warnings,
                "errors": read.parsing_result.get("errors", []),
            })
        response = HttpResponse(read_metrics_tsv(reads), content_type="text/tab-separated-values")
        response['Content-Disposition'] = 'attachment; filename="{}-sanger-read-metrics.tsv"'.format(run.id)
        return response
    if kind == "clustal" and run.clustal_file:
        file_path = os.path.join(settings.MEDIA_ROOT, run.clustal_file.name)
        with open(file_path, "rb") as handle:
            response = HttpResponse(handle.read(), content_type="text/plain")
        response['Content-Disposition'] = 'attachment; filename="{}-sanger.clustal"'.format(run.id)
        return response
    raise Http404


@require_member_can_write_or_admin_project_of_plasmid
def sanger_run_decision(request, plasmid_id, run_id):
    if request.method != "POST":
        raise Http404
    try:
        run = SangerVerificationRun.objects.select_related("plasmid").get(id=run_id, plasmid_id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404
    decision = request.POST.get("manual_decision", "")
    if decision not in ("VERIFIED", "REJECTED", "INCONCLUSIVE", "PENDING", ""):
        raise Http404
    effective_date_raw = request.POST.get("manual_decision_effective_date", "")
    if effective_date_raw:
        try:
            effective_date = datetime.strptime(effective_date_raw, "%Y-%m-%d").date()
        except ValueError:
            raise Http404
    else:
        effective_date = timezone.localdate()
    comment = (request.POST.get("manual_decision_comment", "") or "").strip()[:5000]
    with transaction.atomic():
        plasmid = Plasmid.objects.select_for_update().get(id=plasmid_id)
        locked_run = SangerVerificationRun.objects.select_for_update().get(id=run_id, plasmid=plasmid)
        locked_run.manual_decision = decision
        locked_run.manual_decision_by = request.user if decision else None
        locked_run.manual_decision_at = timezone.now() if decision else None
        locked_run.manual_decision_effective_date = effective_date if decision and decision != "PENDING" else None
        locked_run.manual_decision_comment = comment
        locked_run.save()
        if decision == "VERIFIED":
            plasmid.sequencing_state = 2
            plasmid.sequencing_date = effective_date
            if comment:
                plasmid.sequencing_observations = comment[:1000]
            plasmid.save(update_fields=["sequencing_state", "sequencing_date", "sequencing_observations"])
    return redirect("sanger_run_detail", plasmid_id=plasmid_id, run_id=run_id)


@require_member_can_write_or_admin_project_of_plasmid
def sanger_run_delete(request, plasmid_id, run_id):
    if request.method != "POST":
        raise Http404
    try:
        run = SangerVerificationRun.objects.get(id=run_id, plasmid_id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404
    delete_sanger_run_files(run)
    run.delete()
    return redirect("plasmid_align", plasmid_id=plasmid_id)


class PlasmidDelete(DeleteView):
    model = Plasmid

    @method_decorator(require_member_can_write_or_admin_project_of_plasmid)
    def dispatch(self, *args, **kwargs):
        self.extra_context = {
            'user_can_edit_plasmid': True
        }
        return super().dispatch(*args, **kwargs)

    def get_object(self, *args, **kwargs):
        obj = super(PlasmidDelete, self).get_object(*args, **kwargs)
        if not member_can_write_or_admin_plasmid(obj, self.request.user):
            raise PermissionDenied()
        return obj

    def get_success_url(self, **kwargs):
        return reverse('plasmid_deleted')


def plasmid_deleted(request):
    return render(request, 'inventory/plasmid_deleted.html')


def plasmid_validation_initial_from_payload(validation_payload):
    pattern = re.compile(
        r'^(?P<weaver_id>\d+)_(?P<colony_number>\d+)_'
        r'(?P<date>\d{4}-\d{2}-\d{2})_(?P<method>pcr|digest)$',
        re.IGNORECASE
    )
    match = pattern.match(validation_payload)
    if not match:
        return None, "Invalid validation link format"

    try:
        parsed_date = datetime.strptime(match.group('date'), "%Y-%m-%d").date()
    except ValueError:
        return None, "Invalid date in validation link"

    method = match.group('method').lower()
    initial = {
        'weaver_id': int(match.group('weaver_id')),
        'working_colony': int(match.group('colony_number')),
        'ligation_state': 1,
    }

    if method == 'pcr':
        initial.update({
            'method': 'pcr',
            'colonypcr_state': 2,
            'colonypcr_date': parsed_date,
        })
    elif method == 'digest':
        initial.update({
            'method': 'digest',
            'digestion_state': 2,
            'digestion_date': parsed_date,
        })

    return initial, None


@require_current_project_set
def PlasmidValidationFromLink(request, validation_payload):
    initial, error = plasmid_validation_initial_from_payload(validation_payload)
    if error:
        return render(request, 'inventory/general_error.html', {'error': error})

    try:
        plasmid_to_validate = Plasmid.objects.get(idx=initial['weaver_id'])
    except ObjectDoesNotExist:
        raise Http404

    if not member_can_write_or_admin_plasmid(plasmid_to_validate, request.user):
        raise PermissionDenied

    if request.method == 'POST':
        post_data = request.POST.copy()
        post_data['ligation_state'] = str(initial['ligation_state'])
        post_data['working_colony'] = str(initial['working_colony'])
        if initial['method'] == 'pcr':
            post_data['colonypcr_state'] = str(initial['colonypcr_state'])
            post_data['colonypcr_date'] = initial['colonypcr_date'].isoformat()
        elif initial['method'] == 'digest':
            post_data['digestion_state'] = str(initial['digestion_state'])
            post_data['digestion_date'] = initial['digestion_date'].isoformat()

        form = PlasmidValidationForm(post_data or None, request.FILES, instance=plasmid_to_validate)
        if form.is_valid():
            form.save()
            return HttpResponseRedirect(reverse('plasmid', args=(
                plasmid_to_validate.id,)) + '?form_result_plasmidvalidation_edit_success=true')
    else:
        form = PlasmidValidationForm(instance=plasmid_to_validate, initial=initial)

    return render(request, 'inventory/plasmidvalidation_update_form.html',
                  {'form': form, 'plasmid': plasmid_to_validate,
                   'validation_link_payload': validation_payload,
                   'validation_link_method': initial['method'],
                   'user_can_edit_plasmid': member_can_write_or_admin_plasmid(plasmid_to_validate, request.user)})


@require_current_project_set
@require_member_can_write_or_admin_project_of_plasmid
def PlasmidValidationEdit(request, plasmid_id):
    try:
        plasmid_to_validate = Plasmid.objects.get(id=plasmid_id)
    except ObjectDoesNotExist:
        raise Http404

    if request.method == 'POST':
        form = PlasmidValidationForm(request.POST or None, request.FILES, instance=plasmid_to_validate)
        if form.is_valid():
            form.save()
            return HttpResponseRedirect(reverse('plasmid', args=(
                plasmid_to_validate.id,)) + '?form_result_plasmidvalidation_edit_success=true')
    else:
        form = PlasmidValidationForm(instance=plasmid_to_validate)

    return render(request, 'inventory/plasmidvalidation_update_form.html',
                  {'form': form, 'plasmid': plasmid_to_validate,
        'user_can_edit_plasmid': member_can_write_or_admin_plasmid(plasmid_to_validate, request.user)})


@require_current_project_set
def PlasmidValidations(request):
    # process post
    plasmid_massive_action_results = []
    try:
        if request.method == 'POST':
            if 'massive_action_form_submit' in request.POST and 'massive_action_form_action' in request.POST:
                for key, value in request.POST.items():
                    if key.startswith("pidx-"):
                        plasmid = Plasmid.objects.get(idx=int(key.split("-")[1]))
                        if member_can_write_or_admin_plasmid(plasmid, request.user):
                            action = request.POST.get('massive_action_form_action')
                            action_id = int(action.split("-")[1])
                            if action.startswith("ligation_state"):
                                plasmid.ligation_state = action_id
                                plasmid_massive_action_results.append(
                                    [plasmid, "Set to " + str(LIGATION_STATES[action_id][1])])
                            if action.startswith("colony_pcr"):
                                # set correct & now
                                plasmid.colonypcr_state = action_id
                                if action_id == 2: # correct
                                    plasmid.colonypcr_date = datetime.now()
                                else:
                                    plasmid.colonypcr_date = None
                                plasmid_massive_action_results.append([plasmid, "Set to Colony PCR - " + str(CHECK_STATES[action_id][1])])
                            elif action.startswith("digestion"):
                                # set correct & now
                                plasmid.digestion_state = action_id
                                if action_id == 2: # correct
                                    plasmid.digestion_date = datetime.now()
                                else:
                                    plasmid.digestion_date = None
                                plasmid_massive_action_results.append([plasmid, "Set to Digestion - " + str(CHECK_STATES[action_id][1])])
                            elif action.startswith("sequencing"):
                                # set correct & now
                                plasmid.sequencing_state = action_id
                                if action_id == 2: # correct
                                    plasmid.sequencing_date = datetime.now()
                                else:
                                    plasmid.sequencing_date = None
                                plasmid_massive_action_results.append([plasmid, "Set to Sequencing - " + str(CHECK_STATES[action_id][1])])
                            plasmid.save()
    except Exception:
        pass

    show_from_all_projects = get_show_from_all_projects(request)
    if show_from_all_projects:
        all_plasmids = Plasmid.objects.filter(project_id__in=get_projects_where_member_can_any(request.user))
    else:
        all_plasmids = Plasmid.objects.filter(project_id=get_current_project_id(request))

    plasmidsToStock = []
    for plasmid in all_plasmids:
        if plasmid.ligation_state == 1 and plasmid.reference_sequence is False and plasmid.colonypcr_state != 1 and plasmid.digestion_state != 1 and plasmid.sequencing_state != 1:
            primary_gs = False
            for gs in plasmid.glycerolstock_set.all():
                if gs.strain.for_primary_gs:
                    primary_gs = True
                    break
            if not primary_gs:
                plasmidsToStock.append(plasmid)

    all_plasmids_not_ref = all_plasmids.exclude(reference_sequence=True)
    all_plasmids_ligated = all_plasmids_not_ref.filter(ligation_state=1)

    context = {
        'show_from_all_projects': show_from_all_projects,
        'CHECK_STATES': dict(CHECK_STATES),
        'LIGATION_STATES': dict(LIGATION_STATES),
        'plasmid_massive_action_results': plasmid_massive_action_results,
        'lists': {
            'waiting': {
                'name': 'Waiting parts or supplies',
                'empty_text': 'waiting parts or supplies',
                'data': all_plasmids_not_ref.filter(ligation_state=0)
            },
            'to_ligate': {
                'name': 'Ligation pending',
                'empty_text': 'to ligate',
                'data': all_plasmids_not_ref.exclude(ligation_state=0).exclude(ligation_state=1)
            },
            'to_colonypcr': {
                'name': 'Colony PCR pending',
                'empty_text': 'to colony PCR',
                'data': all_plasmids_ligated.filter(colonypcr_state=1)
            },
            'to_digest': {
                'name': 'Digestion pending',
                'empty_text': 'to digest',
                'data': all_plasmids_ligated.filter(digestion_state=1)
            },
            'to_sequence': {
                'name': 'Sequencing pending',
                'empty_text': 'to sequence',
                'data': all_plasmids_ligated.filter(sequencing_state=1)
            },
            'to_stock': {
                'name': 'To Stock',
                'empty_text': 'to stock',
                'data': plasmidsToStock
            },
            'reference': {
                'name': 'Reference plasmids',
                'empty_text': 'reference',
                'data': all_plasmids.filter(reference_sequence=True)
            },
        }
    }
    return render(request, 'inventory/plasmidvalidations.html', context)


def plasmid_update_computed_size(plasmid_to_update):
    sequence = grab_seq(plasmid_to_update)

    if sequence[0]:
        if not plasmid_to_update.level is None:
            if plasmid_to_update.level % 2:
                re = RestrictionEnzyme.objects.filter(name="BsmBI")[0]
            else:
                re = RestrictionEnzyme.objects.filter(name="BsaI")[0]

            hits = re_find_cut_positions(sequence[1], re, True, True)

            if len(hits) == 2 and hits[1] > hits[0]:
                plasmid_to_update.insert_computed_size = hits[1] - hits[0] + 4

        plasmid_to_update.computed_size = len(sequence[1])
        plasmid_to_update.save()
        return True
    else:
        plasmid_to_update.insert_computed_size = None
        plasmid_to_update.computed_size = None
        plasmid_to_update.save()
        return False


def grab_seq(plasmid_to_grab_from):
    result, gb_record = seqio_get(plasmid_to_grab_from)
    if result:
        return True, gb_record.seq
    return result, gb_record


def grab_features(plasmid_to_grab_from):
    result, gb_record = seqio_get(plasmid_to_grab_from)
    if result:
        return True, gb_record.features
    return result, gb_record


def grab_features_json(plasmid_to_grab_from):
    result, gb_features = grab_features(plasmid_to_grab_from)
    features = []
    if result:
        for gb_feature in gb_features:
            start = int(gb_feature.location.start)
            end = int(gb_feature.location.end)
            if type(gb_feature.location) is Bio.SeqFeature.CompoundLocation:
                if len(gb_feature.location.parts):
                    # asume partes contiguas
                    start = gb_feature.location.parts[1].start
                    end = gb_feature.location.parts[0].end
            forward = False
            if gb_feature.location.strand:
                forward = True
            features.append({
                'name': gb_feature.qualifiers['label'][0],
                'type': gb_feature.type,
                'start': start,
                'end': end,
                'forward': forward,
            })
    return json.dumps(features)


def seqio_get(plasmid_to_grab_from):
    name, extension = os.path.splitext(plasmid_to_grab_from.sequence.name)
    format_name = ''
    if extension == '.gb' or extension == '.gbk':
        format_name = "genbank"
    if extension == '.fasta':
        format_name = "fasta"
    if format_name:
        try:
            for gb_record in SeqIO.parse(plasmid_to_grab_from.sequence.path, format_name):
                return True, gb_record
        except ValueError as e:
            try:
                # Create temp file
                fh, abs_path = mkstemp()
                file_path = plasmid_to_grab_from.sequence.path
                with fdopen(fh, 'w') as new_file:
                    with open(file_path) as old_file:
                        for line in old_file:
                            new_line = line
                            if line.startswith("LOCUS"):
                                line_split = []
                                for idx, val in enumerate(line.split()):
                                    if idx == 3:
                                        continue
                                    if idx == 2:
                                        line_split.append(val + " bp")
                                    else:
                                        line_split.append(val)
                                spaces = [12, 13, 11, 16, 10, 12]
                                new_line = ""
                                for idx, val in enumerate(line_split):
                                    if idx >= len(spaces):
                                        new_line = new_line + val + " "
                                    else:
                                        if spaces[idx] > len(val):
                                            new_line = new_line + val + " " * (spaces[idx] - len(val))
                                        else:
                                            new_line = new_line + val[:spaces[idx] - 1] + " "
                                new_line = new_line + "\n"
                            new_file.write(line.replace(line, new_line))
                # Copy the file permissions from the old file to the new file
                copymode(file_path, abs_path)
                # Remove original file
                remove(file_path)
                # Move new file
                move(abs_path, file_path)
                # Ready
                for gb_record in SeqIO.parse(file_path, format_name):
                    return True, gb_record
            except ValueError as e:
                return False, 'File bad format: ' + e.__str__()
            except FileNotFoundError as e:
                return False, 'File not found: ' + e.__str__()
        except AttributeError as e:
            return False, e.__str__()
    return False, 'Unsupported file extension'


def re_digestion_fragments(sequence, the_res, is_circular):
    ordered_results = []
    for the_re in the_res:
        cut_positions = re_find_cut_positions(sequence, the_re, is_circular, True)
        for cp in cut_positions:
            ordered_results.append((the_re, cp))
    ordered_results.sort(key=lambda tup: tup[1])

    fragments = []
    prev_fragment = None
    last_element = ordered_results[len(ordered_results) - 1]
    for ordered_result in ordered_results:
        if len(fragments) == 0:
            # first item
            if is_circular:
                prev_fragment = {
                    'end': last_element[1],
                    'right': last_element[0]
                }
            else:
                prev_fragment = {
                    'end': 0,
                    'right': 'None'
                }
        fragment = {
            'start': prev_fragment['end'],
            'end': ordered_result[1],
            'left': prev_fragment['right'],
            'right': ordered_result[0],
            'length': ordered_result[1] - prev_fragment['end'],
        }
        if len(fragments) == 0 and is_circular:
            fragment['length'] = ordered_result[1] + len(sequence) - last_element[1]
        fragments.append(fragment)
        prev_fragment = fragment

    fragments.sort(key=lambda dic: dic['length'])
    return fragments


def re_find_cut_fragments(sequence, the_re, is_circular):
    cut_positions = re_find_cut_positions(sequence, the_re, is_circular, True)
    fragments = []
    if cut_positions:
        prev_cut_pos = 0
        for cp in cut_positions:
            fragments.append(cp - prev_cut_pos)
            prev_cut_pos = cp

        last_frag = len(sequence) - cut_positions[len(cut_positions) - 1]
        if is_circular:
            fragments[0] += last_frag
        else:
            fragments.append(last_frag)
    return cut_positions, sorted(fragments)


def re_find_cut_positions(sequence, the_re, is_circular, sort):
    search_results = RestrictionBatch([the_re.name]).search(Seq(sequence), linear=not is_circular)
    found_hits = []
    for key in search_results:
        if str(key) == the_re.name:
            found_hits = search_results[key]
    if sort:
        found_hits = sorted(found_hits)
    return found_hits


def primer_numeric_id(primer):
    return display_primer_id(primer)


def primer_display_name(primer):
    return display_primer_name(primer)


@require_member_can_read_project_of_primer
def primer(request, primer_id):
    try:
        primer_to_detail = Primer.objects.get(id=primer_id)
    except ObjectDoesNotExist:
        raise Http404
    primer_to_detail.display_idx = primer_numeric_id(primer_to_detail)
    primer_to_detail.display_name = primer_display_name(primer_to_detail)
    context = {
        'primer': primer_to_detail,
        'user_can_edit_primer': member_can_write_or_admin_primer(primer_to_detail, request.user)
    }
    return render(request, 'inventory/primer.html', context)


def primers(request):
    show_from_all_projects = get_show_from_all_projects(request)
    if show_from_all_projects:
        primers = visible_primers_for_user(request.user)
    else:
        primers = visible_primers_for_user(request.user)

    primers = sorted(primers, key=lambda primer: (
        primer_numeric_id(primer) if primer_numeric_id(primer) is not None else 10**9,
        primer.name or ""
    ))
    for primer in primers:
        primer.display_idx = primer_numeric_id(primer)
        primer.display_name = primer_display_name(primer)
        primer.can_edit = member_can_write_or_admin_primer(primer, request.user)
    context = {
        'primers': primers,
        'show_from_all_projects': show_from_all_projects,
        'on_current_project_member_can_write_or_admin': on_current_project_member_can_write_or_admin(request)
    }
    return render(request, 'inventory/primers.html', context)


@require_current_project_set
@require_member_can_write_or_admin_current_project
def primer_import(request):
    result = None
    form = PrimerBatchUploadForm(request.POST or None, request.FILES or None)
    current_project = get_current_project(request)

    if request.method == "POST" and form.is_valid():
        try:
            fasta_text = request.FILES["fasta_file"].read().decode("utf-8-sig")
            result = import_primers_from_fasta(
                StringIO(fasta_text),
                current_project,
                update_existing=form.cleaned_data["update_existing"],
                name_source=form.cleaned_data["name_source"],
                default_direction=form.cleaned_data["default_direction"],
            )
        except UnicodeDecodeError:
            form.add_error("fasta_file", "Could not read this file as UTF-8 text.")
        except PrimerImportError as error:
            form.add_error("fasta_file", str(error))

    context = {
        "form": form,
        "result": result,
        "current_project": current_project,
    }
    return render(request, "inventory/primer_import.html", context)


class PrimerCreate(CreateView):
    model = Primer
    fields = '__all__'
    template_name_suffix = '_create_form'

    @method_decorator(require_member_can_write_or_admin_current_project)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_success_url(self, **kwargs):
        return reverse('primer', args=(self.object.id,)) + '?form_result_primer_create_success=true'


class PrimerEdit(UpdateView):
    model = Primer
    fields = '__all__'
    template_name_suffix = '_update_form'

    @method_decorator(require_member_can_write_or_admin_project_of_primer)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_success_url(self, **kwargs):
        return reverse('primer', args=(self.object.id,)) + '?form_result_primer_edit_success=true'


class PrimerDelete(DeleteView):
    model = Primer

    @method_decorator(require_member_can_write_or_admin_project_of_primer)
    def dispatch(self, *args, **kwargs):
        return super().dispatch(*args, **kwargs)

    def get_success_url(self, **kwargs):
        return reverse('primers') + '?form_result_object_deleted=true'


@require_member_can_read_project_of_primer
def primer_label(request, primer_id):
    try:
        primer_to_label = Primer.objects.get(id=primer_id)
    except ObjectDoesNotExist:
        raise Http404

    context = {
        'primer': primer_to_label,
        'date': datetime.now().date(),
    }
    return render(request, 'inventory/primer_label.html', context)


def ServicesStats(request):
    if Stats.objects.all():
        stats = Stats.objects.all()[0]
    else:
        stats = Stats()
        stats.save()

    context = {
        'error': 'No data'
    }

    if request.method == 'POST' and 'refresh_stats' in request.POST:
        plasmids_by_month = {'date': [], 'plasmid_month_count': []}
        current_year = ''
        current_month = ''
        current_month_count = 0
        plasmid_count = 0
        plasmids_with_sequence = 0
        plasmids_with_gs = 0
        plasmids_by_type = {}
        plasmids_by_level = {}
        plasmids_ordered = Plasmid.objects.order_by('created_on')
        for plasmid in plasmids_ordered:
            plasmid_count += 1
            year = plasmid.created_on.year
            month = plasmid.created_on.month
            if current_year == '' or current_month == '':
                current_year = year
                current_month = month
            if current_year == year and current_month == month:
                current_month_count += 1
            else:
                # save
                plasmids_by_month['date'].append(json_serial(datetime(current_year, current_month, 1)))
                plasmids_by_month['plasmid_month_count'].append(current_month_count)
                # update current
                current_year = year
                current_month = month
                current_month_count = 1
            if plasmid.computed_size:
                plasmids_with_sequence += 1
            if plasmid.glycerolstock_set.all().count():
                plasmids_with_gs += 1

            key = str(plasmid.type)
            if key not in plasmids_by_type:
                plasmids_by_type[key] = 0
            plasmids_by_type[key] += 1

            key = "Level " + str(plasmid.level)
            if key not in plasmids_by_level:
                plasmids_by_level[key] = 0
            plasmids_by_level[key] += 1

        # save last
        plasmids_by_month['date'].append(json_serial(datetime(current_year, current_month, 1)))
        plasmids_by_month['plasmid_month_count'].append(current_month_count)

        stats.plasmids_by_month = plasmids_by_month
        stats.plasmids_with_sequence = {'values': [plasmids_with_sequence, plasmid_count - plasmids_with_sequence],
                                        'names': ['With sequence', 'Without sequence']}
        stats.plasmids_with_gs = {'values': [plasmids_with_gs, plasmid_count - plasmids_with_gs],
                                  'names': ['With GStock', 'Without Gstock']}

        stats.plasmids_by_type = {'values': list(plasmids_by_type.values()), 'names': list(plasmids_by_type.keys())}
        stats.plasmids_by_level = {'values': list(plasmids_by_level.values()), 'names': list(plasmids_by_level.keys())}
        stats.last_update = date.today()
        stats.plasmid_count = plasmid_count
        stats.save()

    if stats.plasmids_by_month:
        df = pd.DataFrame(data=stats.plasmids_by_month)
        fig_plasmid_month_count = px.line(df, x="date", y="plasmid_month_count", text="plasmid_month_count",
                                          title="Plasmid creation",
                                          labels={'date': 'Date', 'plasmid_month_count': '# created plasmids'})
        fig_plasmid_month_count.update_traces(textposition="top center")

        df = pd.DataFrame(data=stats.plasmids_with_sequence)
        fig_plasmids_with_sequence = px.pie(df, values='values', names='names', hole=.3, title="Plasmids Sequence")

        df = pd.DataFrame(data=stats.plasmids_with_gs)
        fig_plasmids_with_gs = px.pie(df, values='values', names='names', hole=.3, title="Plasmids GStock")

        df = pd.DataFrame(data=stats.plasmids_by_type)
        fig_plasmids_by_type = px.pie(df, values='values', names='names', hole=.3, title="Type")

        df = pd.DataFrame(data=stats.plasmids_by_level)
        fig_plasmids_by_level = px.pie(df, values='values', names='names', hole=.3, title="Level")

        context = {
            'fig_plasmid_month_count': fig_plasmid_month_count.to_html(),
            'fig_plasmids_with_sequence': fig_plasmids_with_sequence.to_html(),
            'fig_plasmids_with_gs': fig_plasmids_with_gs.to_html(),
            'fig_plasmids_by_type': fig_plasmids_by_type.to_html(),
            'fig_plasmids_by_level': fig_plasmids_by_level.to_html(),
            'last_update': stats.last_update,
            'plasmid_count': stats.plasmid_count
        }
    return render(request, 'inventory/services/stats/stats.html', context)


def ServicesGtr(request):
    return render(request, 'inventory/services/gtr/gtr.html')


def ServicesL0d(request):
    context={
        'csrf_token': django.middleware.csrf.get_token(request),
    }
    return render(request, 'inventory/services/l0d/l0d.html', context)


def batch_print_tokens(raw_identifiers):
    return [token.strip() for token in re.split(r"[\n,;]+", raw_identifiers or "") if token.strip()]


def looks_like_uuid(value):
    try:
        uuid.UUID(str(value))
        return True
    except ValueError:
        return False


def find_batch_plasmids(tokens, user):
    queryset = Plasmid.objects.filter(project__in=get_projects_where_member_can_any(user))
    objects = []
    missing = []
    for token in tokens:
        query = Q(name=token) | Q(qr_id=token)
        if token.isdigit():
            query |= Q(idx=int(token))
        if looks_like_uuid(token):
            query |= Q(id=token)
        match = queryset.filter(query).first()
        if match:
            objects.append(match)
        else:
            missing.append(token)
    return objects, missing


def find_batch_glycerolstocks(tokens, user):
    queryset = GlycerolStock.objects.filter(project__in=get_projects_where_member_can_any(user)).select_related(
        "plasmid", "strain", "box", "box__location"
    )
    objects = []
    missing = []
    for token in tokens:
        query = Q(qr_id=token)
        if looks_like_uuid(token):
            query |= Q(id=token)
        if token.isdigit():
            query |= Q(plasmid__idx=int(token))
        match = queryset.filter(query).first()
        if match:
            if match.plasmid:
                match.resistantes_human = resistantes_human(match.plasmid.selectable_markers, True)
            else:
                match.resistantes_human = "None"
            match.resistantes_strain_human = resistantes_human(match.strain.selectable_markers, True)
            objects.append(match)
        else:
            missing.append(token)
    return objects, missing


def find_batch_label(label_type, token, user):
    if label_type == "plasmids":
        matches, missing = find_batch_plasmids([token], user)
    elif label_type == "glycerolstocks":
        matches, missing = find_batch_glycerolstocks([token], user)
    else:
        return None, "Unknown label type"
    if matches:
        return matches[0], ""
    return None, missing[0] if missing else token


def parse_batch_print_date(raw_value):
    try:
        return date.fromisoformat(raw_value)
    except (TypeError, ValueError):
        return date.today()


def parse_batch_print_concentration(raw_value):
    try:
        value = float(raw_value)
    except (TypeError, ValueError):
        return ""
    if value <= 0:
        return ""
    return value


def parse_batch_print_colony(raw_value):
    colony = (raw_value or "").strip()
    if colony.lower().startswith("c"):
        colony = colony[1:].strip()
    return colony


def batch_print_rows_from_post(post_data):
    label_types = post_data.getlist("label_type")
    identifiers = post_data.getlist("identifier")
    colonies = post_data.getlist("colony")
    dates = post_data.getlist("date")
    concentrations = post_data.getlist("concentration")
    rows = []
    total = max(len(label_types), len(identifiers), len(colonies), len(dates), len(concentrations))
    for index in range(total):
        row = {
            "label_type": label_types[index] if index < len(label_types) else "plasmids",
            "identifier": identifiers[index].strip() if index < len(identifiers) else "",
            "colony": parse_batch_print_colony(colonies[index]) if index < len(colonies) else "",
            "date": dates[index] if index < len(dates) else date.today().isoformat(),
            "concentration": concentrations[index] if index < len(concentrations) else "",
        }
        if row["identifier"]:
            rows.append(row)
    return rows


def ServicesBatchPrints(request):
    labels = []
    missing = []
    searched = False
    rows = [{
        "label_type": "plasmids",
        "identifier": "",
        "colony": "",
        "date": date.today().isoformat(),
        "concentration": "",
    }]
    if request.method == "POST":
        searched = True
        rows = batch_print_rows_from_post(request.POST)
        for row in rows:
            obj, missing_token = find_batch_label(row["label_type"], row["identifier"], request.user)
            if obj:
                labels.append({
                    "kind": row["label_type"],
                    "object": obj,
                    "colony": row["colony"],
                    "date": parse_batch_print_date(row["date"]),
                    "concentration": parse_batch_print_concentration(row["concentration"]),
                })
            else:
                missing.append("{} ({})".format(row["identifier"], missing_token))
        if not rows:
            rows = [{
                "label_type": "plasmids",
                "identifier": "",
                "colony": "",
                "date": date.today().isoformat(),
                "concentration": "",
            }]

    context = {
        "rows": rows,
        "label_types": BatchPrintsForm.LABEL_TYPES,
        "labels": labels,
        "missing": missing,
        "searched": searched,
        "today": date.today().isoformat(),
    }
    return render(request, "inventory/services/batch_prints/batch_prints.html", context)


def ServicesPcr(request):
    form = ServicesPCRForm(request.POST or None, user=request.user)
    results = []
    skipped = []
    primer_f = None
    primer_r = None
    primer_complementarity = None

    if request.method == "POST" and form.is_valid():
        primer_f = form.cleaned_data["primer_f"]
        primer_r = form.cleaned_data["primer_r"]
        primer_complementarity = primer_pair_complementarity(primer_f, primer_r)
        min_product_size = form.cleaned_data["min_product_size"]
        max_product_size = form.cleaned_data["max_product_size"]
        plasmids = Plasmid.objects.filter(
            project__in=get_projects_where_member_can_any(request.user)
        ).order_by("name")

        for plasmid_to_scan in plasmids:
            sequence = grab_seq(plasmid_to_scan)
            if not sequence[0]:
                skipped.append({
                    "plasmid": plasmid_to_scan,
                    "reason": sequence[1],
                })
                continue

            for amplicon in primer_pair_amplicons(
                    str(sequence[1]),
                    primer_f,
                    primer_r,
                    min_product_size=min_product_size,
                    max_product_size=max_product_size,
            ):
                results.append({
                    "plasmid": plasmid_to_scan,
                    "project": plasmid_to_scan.project,
                    "amplicon": amplicon,
                })

        results.sort(key=lambda result: (
            result["amplicon"]["product_size"],
            result["plasmid"].name or "",
        ))

    context = {
        "form": form,
        "results": results,
        "skipped": skipped,
        "searched": request.method == "POST" and form.is_valid(),
        "primer_f": primer_f,
        "primer_r": primer_r,
        "primer_complementarity": primer_complementarity,
    }
    return render(request, 'inventory/services/pcr/pcr.html', context)


def run_local_blast(request, context, record, project_id='a', short_blast=False):
    visible_projects = get_projects_where_member_can_any(request.user)
    if project_id == 'a':
        plasmids = Plasmid.objects.filter(project__in=visible_projects)
    else:
        plasmids = Plasmid.objects.filter(project__in=visible_projects, project=project_id)

    subjects = []
    context['not_considered_subjects'] = []
    for plasmid in plasmids:
        try:
            seqio_get_result = seqio_get(plasmid)
            if seqio_get_result[0]:
                seqio_get_result[1].id = plasmid.id
                seqio_get_result[1].name = plasmid.name
                subjects.append(make_circular([seqio_get_result[1]])[0])
            else:
                context['not_considered_subjects'].append((plasmid, 'No sequence file'))
        except Exception as e:
            context['not_considered_subjects'].append((plasmid, e))

    if record and subjects:
        context['query'] = record
        context['short_blast'] = "Yes" if short_blast else "No"
        queries = make_linear([record])
        blast = BioBlast(subjects, queries)
        context['results'] = run_pyblast_compat(
            lambda: blast.blastn_short() if short_blast else blast.blastn()
        )
        for result in context['results']:
            result['alignment'] = Bio.Align.MultipleSeqAlignment([
                SeqRecord(Seq(result['meta']['query seq']), id=result['query']['name']),
                SeqRecord(Seq(result['meta']['subject seq']), id=result['subject']['name']),
            ]).__format__('clustal')
    elif not context.get('error'):
        context['error'] = 'Error while parsing input sequence'


def fasta_record_from_text(text, name='Query'):
    sequence_text = str(text or '').strip()
    if not sequence_text:
        return None
    if ">" not in sequence_text:
        sequence_text = f">{name}\n{sequence_text}\n"
    records = fasta_records_from_text(sequence_text)
    return records[0] if records else None


def ServicesBlast(request):
    context = {}
    project_choices = [('a', 'All')]
    for project in get_projects_where_member_can_any(request.user):
        project_choices.append((project.id, project.name),)

    if request.method == 'GET' and request.GET.get('sequence'):
        record = fasta_record_from_text(request.GET.get('sequence'), request.GET.get('name') or 'Amplicon')
        if record is None:
            context['error'] = 'Input sequence not in FASTA format'
        run_local_blast(
            request,
            context,
            record,
            project_id=request.GET.get('project', 'a'),
            short_blast=request.GET.get('short_blast', '').lower() in ('1', 'true', 'yes'),
        )
    elif request.method == 'POST':
        form = BlastSequenceInput(project_choices, request.POST, request.FILES)
        if form.is_valid():
            if request.POST.get('fasta_sequence'):
                record = fasta_record_from_text(request.POST.get('fasta_sequence'))
                if record is None:
                    context['error'] = 'Input sequence not in FASTA format'
            else:
                if request.FILES["fasta_file"]:
                    fasta_file = request.FILES["fasta_file"]
                    f = tempfile.NamedTemporaryFile(delete=False)
                    for chunk in fasta_file.chunks():
                        f.write(chunk)
                    f.close()
                    record = SeqIO.read(f.name, "fasta")
                else:
                    record = None
                    context['error'] = 'No input sequence'

            run_local_blast(
                request,
                context,
                record,
                project_id=request.POST.get('project'),
                short_blast=bool(request.POST.get('short_blast')),
            )
    else:
        context['form'] = BlastSequenceInput(project_choices)

    return render(request, 'inventory/services/blast/blast.html', context)


def get_plasmid_type_id(plasmid):
    if plasmid.type:
        return plasmid.type.id
    return None


def api_plasmid_get_fasta(plasmid):
    if plasmid:
        sequence = ""
        result = grab_seq(plasmid)
        if result[0]:
            sequence = result[1]
        return JsonResponse({
            'name': str(plasmid),
            'id': str(plasmid.id),
            'idx': str(plasmid.idx),
            'seq': str(sequence)
        })
    return JsonResponse({
        'error': 'Plasmid not found'
    })


def api_plasmid_get_fasta_by_idx(request, idx):
    return api_plasmid_get_fasta(Plasmid.objects.get(idx=idx))


def api_plasmid_get_fasta_by_name(request, name):
    return api_plasmid_get_fasta(Plasmid.objects.get(name=name))


@require_member_can_any_current_project
def api_plasmids(request):
    output = []
    level_from_table_filters = 0
    level_to_table_filters = 0
    if get_show_from_all_projects(request):
        plasmids = Plasmid.objects.filter(project_id__in=get_projects_where_member_can_any(request.user)).order_by(
            'name')
    else:
        plasmids = Plasmid.objects.filter(project_id=get_current_project_id(request)).order_by('name')
    for plasmid in plasmids:
        output.append({
            'cn': plasmid.__str__(),
            'n': plasmid.name,
            'l': plasmid.level,
            't': get_plasmid_type_id(plasmid),
            'i': plasmid.id,
            'ix': str(plasmid.idx),
            'hs': bool(plasmid.sequence),
            'c': plasmid.computed_size,
            'ic': plasmid.insert_computed_size,
            'cs': plasmid.get_check_state(),
            'r': plasmid.recommended_enzyme_for_create(),
            'sm': " + ".join(list(plasmid.selectable_markers.all().values_list('three_letter_code', flat=True))),
            'p': member_can_write_or_admin_plasmid(plasmid, request.user),
            'wc': plasmid.working_colony_text(),
            'lc': plasmid.ligation_concentration(),
            'd': plasmid.description,
            'iu': plasmid.intended_use
        })
        if plasmid.level:
            if plasmid.level > level_to_table_filters:
                level_to_table_filters = plasmid.level
            if plasmid.level < level_from_table_filters:
                level_from_table_filters = plasmid.level
    context = {
        'table_filters': get_table_filters(level_from_table_filters, level_to_table_filters),
        'plasmids': output,
        'csrf_token': django.middleware.csrf.get_token(request),
        'RESTRICTION_ENZYMES': list(RestrictionEnzyme.objects.values()),
    }
    return JsonResponse(context, safe=False)


@require_member_can_any_current_project
def api_glycerolstocks(request):
    output = []
    level_from_table_filters = 0
    level_to_table_filters = 0
    if get_show_from_all_projects(request):
        glycerolstocks = GlycerolStock.objects.filter(
            project_id__in=get_projects_where_member_can_any(request.user)).order_by('strain', 'plasmid')
    else:
        glycerolstocks = GlycerolStock.objects.filter(project_id=get_current_project_id(request)).order_by('strain',
                                                                                                           'plasmid')
    for glycerolstock in glycerolstocks:
        pi = ""
        pix = ""
        pn = ""
        pt = ""
        pl = ""
        pcs = ''
        if glycerolstock.plasmid:
            if glycerolstock.plasmid.level:
                if glycerolstock.plasmid.level > level_to_table_filters:
                    level_to_table_filters = glycerolstock.plasmid.level
                if glycerolstock.plasmid.level < level_from_table_filters:
                    level_from_table_filters = glycerolstock.plasmid.level
            pi = glycerolstock.plasmid.id
            pix = glycerolstock.plasmid.idx
            pn = glycerolstock.plasmid.name
            if glycerolstock.plasmid.type:
                pt = glycerolstock.plasmid.type.id
            pl = glycerolstock.plasmid.level
            pcs = glycerolstock.plasmid.get_check_state()
        bn = ""
        bl = ""
        if glycerolstock.box:
            bn = glycerolstock.box.name
            bl = str(glycerolstock.box.location)
        output.append({
            'i': glycerolstock.id,
            'pi': pi,
            'pix': pix,
            'pcs': pcs,
            'pn': pn,
            'pt': pt,
            'pl': pl,
            's': str(glycerolstock.strain),
            'bc': glycerolstock.box_column,
            'br': glycerolstock.box_row,
            'bn': bn,
            'bl': bl,
            'p': member_can_write_or_admin_gs(glycerolstock, request.user)
        })
    context = {
        'table_filters': get_table_filters(level_from_table_filters, level_to_table_filters),
        'glycerolstocks': output,
    }
    return JsonResponse(context, safe=False)


def experiments(request):
    projects_with_experiments = []

    for project in get_projects_where_member_can_any(request.user):
        if project.experiment_set.all():
            projects_with_experiments.append(project)

    context = {
        'projects': projects_with_experiments
    }
    return render(request, 'inventory/experiments.html', context)


def _experiment_plasmid_status(plasmid):
    if plasmid.reference_sequence:
        return 'RS'
    if plasmid.is_validated():
        return 'V'
    if plasmid.ligation_state != 1:
        return 'UC'
    return 'NV'


def _experiment_plasmid_node(plasmid, request):
    plasmid_url = reverse('plasmid', kwargs={'plasmid_id': plasmid.id})
    return {
        'uuid': str(plasmid.id),
        'weaver_id': plasmid.idx,
        'name': plasmid.name,
        'plasmid_type': str(plasmid.type) if plasmid.type else '',
        'type_id': plasmid.type.id if plasmid.type else None,
        'level': plasmid.level,
        'parts': [],
        'parent': [],
        'status': _experiment_plasmid_status(plasmid),
        'colony': plasmid.working_colony,
        'url': request.build_absolute_uri(plasmid_url),
        'ligation_raw': plasmid.ligation_raw(),
        'ready_to_build': False,
    }


def _collect_experiment_plasmid(plasmid, request, nodes, visiting=None):
    if visiting is None:
        visiting = set()
    if plasmid.idx is None or plasmid.idx in visiting:
        return

    visiting.add(plasmid.idx)
    if plasmid.idx not in nodes:
        nodes[plasmid.idx] = _experiment_plasmid_node(plasmid, request)

    child_ids = []
    children = []
    if plasmid.backbone:
        children.append(plasmid.backbone)
    children.extend(list(plasmid.inserts.all()))

    for child in children:
        if child.idx is None:
            continue
        child_ids.append(child.idx)
        _collect_experiment_plasmid(child, request, nodes, visiting.copy())
        if child.idx in nodes and plasmid.idx not in nodes[child.idx]['parent']:
            nodes[child.idx]['parent'].append(plasmid.idx)

    for child_id in child_ids:
        if child_id not in nodes[plasmid.idx]['parts']:
            nodes[plasmid.idx]['parts'].append(child_id)


def _experiment_map_stats(nodes):
    total = len(nodes)
    validated = len([node for node in nodes.values() if node['status'] == 'V'])
    reference = len([node for node in nodes.values() if node['status'] == 'RS'])
    pending_nodes = [
        node for node in nodes.values()
        if node['status'] not in ('V', 'RS') and node['level'] is not None
    ]
    ready_to_build = 0
    blocked = 0

    for node in pending_nodes:
        dependencies_ready = all(
            nodes[part_id]['status'] in ('V', 'RS')
            for part_id in node['parts']
            if part_id in nodes
        )
        node['ready_to_build'] = dependencies_ready
        if dependencies_ready:
            ready_to_build += 1
        else:
            blocked += 1

    return {
        'total': total,
        'validated': validated,
        'reference': reference,
        'pending': len(pending_nodes),
        'ready_to_build': ready_to_build,
        'blocked': blocked,
        'progress': round((validated / total) * 100) if total else 0,
    }


def api_experiments_map(request):
    projects = []
    for project in get_projects_where_member_can_any(request.user):
        experiments_output = []
        for experiment in project.experiment_set.all():
            nodes = {}
            root_ids = []
            for plasmid in experiment.plasmids.all():
                if plasmid.idx is None:
                    continue
                root_ids.append(plasmid.idx)
                _collect_experiment_plasmid(plasmid, request, nodes)

            experiments_output.append({
                'id': experiment.id,
                'name': experiment.name,
                'description': experiment.description,
                'root_ids': root_ids,
                'stats': _experiment_map_stats(nodes),
                'plasmids': list(nodes.values()),
            })

        if experiments_output:
            projects.append({
                'id': project.id,
                'name': str(project),
                'experiments': experiments_output,
            })

    return JsonResponse({'projects': projects})


def createEnzymeFromName(enzyme_name):
    if enzyme_name in rest_dict:
        return RestrictionEnzyme(name=enzyme_name)
    return None


def api_parts(request, enzyme_name, assembly_standard):
    parts = []
    api_error = ""

    the_re = createEnzymeFromName(enzyme_name)
    if not the_re:
        api_error = 'API Error / Restriction enzyme not found (' + enzyme_name + ')'
    else:
        for plasmid in Plasmid.objects.filter(
                project__in=get_projects_where_member_can_any(request.user)
        ).order_by('name'):
            if (plasmid.level is not None and plasmid.type is not None and plasmid.type.id == 1 and
                    (assembly_standard == 'loop' or assembly_standard == 'ytk')
                    and (
                            (enzyme_name == 'BsaI' and plasmid.level % 2 == 0) or
                            (enzyme_name == 'SapI' and plasmid.level % 2 == 1) or
                            (enzyme_name == 'BsmBI' and plasmid.level % 2 == 1)
                    )):
                # make sure use correct enzyme at loop
                continue
            plasmid_grab_seq = grab_seq(plasmid)
            if plasmid.level is not None and plasmid_grab_seq[0]:
                found_cut_positions = re_find_cut_positions(plasmid_grab_seq[1], the_re, True, True)
                if len(found_cut_positions) == 2:
                    length = found_cut_positions[1] - found_cut_positions[0]
                    oh5 = str(plasmid_grab_seq[1])[
                          found_cut_positions[0] - 1:found_cut_positions[0] + abs(the_re.fcut - the_re.rcut) - 1]
                    oh3 = str(plasmid_grab_seq[1])[
                          found_cut_positions[1] - 1:found_cut_positions[1] + abs(the_re.fcut - the_re.rcut) - 1]
                    if get_plasmid_type_id(plasmid) == 1:
                        # receiver
                        ohtmp = oh3
                        oh3 = oh5
                        oh5 = ohtmp
                        length = len(str(plasmid_grab_seq[1])) - length
                    parts.append({
                        'n': plasmid.__str__(),
                        'd': plasmid.description,
                        'l': plasmid.level,
                        't': get_plasmid_type_id(plasmid),
                        'i': plasmid.id,
                        'len': length,
                        'o5': oh5,
                        'o3': oh3,
                    })
    context = {
        'error': api_error,
        'parts': parts,
        'csrf_token': django.middleware.csrf.get_token(request),
    }
    return JsonResponse(context, safe=False)


def api_fidelity_calc(request, enzyme, ohs):
    url = 'https://ligasefidelity.neb.com/viewset/run.cgi'
    context = {
        'error': 'Bad parameters'
    }
    if enzyme == 'sapi':
        page = requests.post(url, {
            'ohlen': 3,
            'dataset': 'b3-SapI-37_16_cycling',
            'olist': ','.join(ohs.split('-'))
        })
        soup = BeautifulSoup(page.content, "html.parser")
        ligation_fidelity = re.findall(r'\d+', soup.find_all('div', class_="estimated-fidelity")[0].text)[0]
        ligation_frequency_matrix = soup.find_all('div', class_="tool-results-header")[0].findNext('pre')
        ligation_frequency_matrix_html = ligation_frequency_matrix.__str__() + ligation_frequency_matrix.find_next_siblings('table')[0].__str__()
        context = {
            'fidelity': ligation_fidelity,
            'ligation_frequency_matrix_html': ligation_frequency_matrix_html
        }
    return JsonResponse(context)
