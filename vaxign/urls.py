from django.urls import path

from vaxign.views import views
from vaxign.views import queries
from vaxign.views import runs
from vaxign.views import projects
from vaxign.views import vaxitop

from . import apis

urlpatterns = [
    path('', views.index, {'name':'index'}),
    path('dynamic', views.index, {'name':'dynamic'}),
    path('precompute', views.index, {'name':'precompute'}),
    
    path('aucs', views.aucs),
    path('disclaimer', views.disclaimer),
    path('docs', views.docs),
    path('faqs', views.faqs),
    path('tutorial', views.tutorial),
    path('updates', views.updates),
    path('stats', views.stats),
    
    path('login', views.login),
    path('logout', views.logout),
    
    path('query/<str:queryID>', queries.results),
    path('query/<str:queryID>/results', queries.results),
    path('query/<str:queryID>/ortholog', queries.ortholog),
    path('query/<str:queryID>/protein/<str:seqID>', queries.protein),
    path('query/<str:queryID>/protein/<str:seqID>/blast/<str:organism>', queries.blast),
    
    path('query/<str:queryID>/vaxitop', vaxitop.index),
    path('query/<str:queryID>/protein/<str:seqID>/vaxitop', vaxitop.protein),
    path('query/<str:queryID>/protein/<str:seqID>/iedb/search', vaxitop.iedb_search),
    
    path('run/<str:queryID>', runs.index),
    path('run/<str:queryID>/vaxitop', runs.vaxitop),
    
    path('vaxign-ml', views.vaxignml),
    
    path('vaxitop',views.vaxitop),
    
    path('project', projects.index),
    path('project/add', projects.add),
    path('project/<str:projectID>', projects.open),
    path('project/<str:projectID>/edit', projects.edit),
    path('project/<str:projectID>/remove', projects.remove),
    path('project/<str:projectID>/curator/<email>/add', projects.curator, {'type':'add'}),
    path('project/<str:projectID>/curator/<email>/remove', projects.curator, {'type':'remove'}),
    path('project/<str:projectID>/query/run', projects.run),
    path('project/<str:projectID>/query/<str:queryID>/setting', projects.querySetting),
    path('project/<str:projectID>/query/<str:queryID>/remove', projects.queryRemove),
    path('project/<str:projectID>/query/<str:queryID>/public', projects.queryPermission, {'type':1}),
    path('project/<str:projectID>/query/<str:queryID>/private', projects.queryPermission, {'type':0}),
]

# API patterns
urlpatterns += [
    path('api/t_vaxign_query/group', apis.t_vaxign_query_group_all),
    path('api/t_vaxign_query/group/<str:genomeGroup>', apis.t_vaxign_query_group),
    path('api/t_vaxign_query/ortholog/exclude/<str:queryID>', apis.t_vaxign_query_ortholog_exclude),
    path('api/t_user_query/all', apis.t_user_query_all),
    
    path('api/t_vaxign_allele_group/<str:species>/I', apis.t_vaxign_allele_group, {'mhc':'I'}),
    path('api/t_vaxign_allele_group/<str:species>/II', apis.t_vaxign_allele_group, {'mhc':'II'}),
    
    path('api/t_vaxign_mast_results/<str:queryID>', apis.t_vaxign_mast_results),
    path('api/t_vaxign_mast_results/<str:queryID>/<str:seqID>', apis.t_vaxign_mast_results_one_seq),
    
    path('api/iedb_epitope_t_cell_linear/<str:seqID>', apis.iedb_epitope_t_cell_linear),
]
