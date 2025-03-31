---
layout: page
title: FE-Project publications
---

## General References
{% for paper in site.data.publications.general %}
1. {{ paper.authors }}, {{ paper.year }}, {{ paper.title }}. {{ paper.journal }}, {{ paper.volume }}, {{ paper.pages }}, [doi:{{ paper.doi }}](https://doi.org/{{ paper.doi }})
{% endfor %}

## Publications

{% assign year = "" %}
{% for paper in site.data.publications.papers %}
{%   unless paper.lang == "ja" %}
{%     unless paper.year == year %}
{%       assign year = paper.year %}
### {{ year }}
{%     endunless %}
1. {{ paper.authors }}, {{ paper.year }}, {{ paper.title }}. {{ paper.journal }},{% if paper.volume %} {{ paper.volume }}, {% endif %} {% if paper.publisher %} {{ paper.publisher }}, {% endif %} {% if paper.pages %} {{ paper.pages }}, {% endif %} {% if paper.doi %} [doi:{{ paper.doi }}](https://doi.org/{{ paper.doi }}){:target="_blank"} {% endif %}
{%   endunless %}
{% endfor %}
