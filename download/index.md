---
layout: page
title: Download
---

* [Source code](#source-code)
* [Documents](#documents)



## Source code

### Latest version

{% assign dlsite = "https://github.com/ywkawai/FE-Project/archive/refs/tags" %}

{% assign ver = site.data.currentversion.version %}
Version : {{ ver }}

* update: {{ site.data.currentversion.update }}
* source file: [FE-Project-{{ ver }}.tar.gz]({{ dlsite }}/FE-Project-{{ ver }}.tar.gz)
  * [GitHub repository](https://github.com/ywkawai/FE-Project/releases/tag/v{{ site.data.currentversion.version }})

### Old versions

{% for ver in site.data.oldversions %}
{% if ver.oldname %}
* [{{ ver.name }}.tar.gz]({{ dlsite }}/{{ ver.name }}.tar.gz) ({{ ver.oldname }}), {{ ver.date }} <br/>
{% else %}
* [{{ ver.name }}.tar.gz]({{ dlsite }}/{{ ver.name }}.tar.gz), {{ ver.date }} <br/>
{% endif %}
{% endfor %}


### History

{% for ver in site.data.history %}
<details>
  <summary>
    {{ ver.ver }} ( {{ ver.type }} from {{ ver.prev }} )
  </summary>
  <div>
    <ul>
    <li>FElib<ul>{% for item in ver.item_felib %}
      <li>{{ item }}</li>
  {% endfor %}</ul></li>
    </ul>
    <ul>
    <li>Sample programs using FElib<ul>{% for item in ver.item_sample %}
      <li>{{ item }}</li>
  {% endfor %}</ul></li>
    </ul>
    <ul>
    <li>Models using FElib<ul>{% for item in ver.item_model %}
      <li>{{ item }}</li>
  {% endfor %}</ul></li>
    </ul>
  </div>
</details>
{% endfor %}


## Documents
The documents for FE-Project are available at <a href="{{ '/documents/' | relative_url }}">FE-Project Document page</a>


