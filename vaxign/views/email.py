from datetime import datetime
from django.core.mail import send_mail

from django.conf import settings

import logging
logger = logging.getLogger('console')


def email_submission(queryID, email_address, email_body):
    to = [email_address]
    subject = str.format("Vaxign job {} submitted", queryID)
    html_message = str.format("""\
<pre>
<strong>Your Vaxign job submitted</strong>

<a href="http://violinet.org/vaxign2/query/{}/results">Click here</a> to browse your job

==================================================
Job ID
==================================================
{}

==================================================
Submission Time
==================================================
{}

{}

</pre>
""", queryID, queryID, datetime.now(), email_body)
    
    send_mail(subject, '', settings.EMAIL_HOST_USER, to, html_message=html_message)
    
