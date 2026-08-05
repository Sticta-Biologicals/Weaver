from inventory.models import Primer
from organization.views import get_projects_where_member_can_any


def visible_primers_for_user(user):
    if not user or not user.is_authenticated:
        return Primer.objects.none()
    return Primer.objects.filter(project__in=get_projects_where_member_can_any(user))

