---
layout: default
title: FE-Project top page
lang: en
---

<!-- Global site tag (gtag.js) - Google Analytics -->
<script async src="https://www.googletagmanager.com/gtag/js?id=G-8KLNNQVBZF"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'G-8KLNNQVBZF');
</script>

# What is FE-Project?

<p>
This project provides a library and sample programs for the discontinuous Galerkin methods (DGM). Futhermore, dynamical cores of atmospheric model with DGM will be provided. 
</p>

![Top Image]({{ '/assets/images/FE-project_top.png' | relative_url }})


# Contributors

- main developers
  - [Yuta Kawai](https://researchmap.jp/ykawai1988/?lang=english)<sup>*1</sup> 
- project design and scientific decision 
  - Yuta Kawai<sup>*1</sup>  and Hirofumi Tomita<sup>*1</sup>  

Affiliation: *1: [RIKEN Center for Computational Science](http://www.r-ccs.riken.jp/en/), Kobe in Japan

# Acknowledgements

This project is supported by 
[the Transformaive Research Areas B: DNA Climate Science](https://dna-climate.org/) (MEXT KAKENHI Grant Number JP20H05731) 
and 
JST AIP Grant Number JPMJCR19U2. 
The model development and numerical experiments are
performed using the Oackbridge-CX supercomputer at the University of Tokyo and the supercomputer Fugaku at RIKEN (Project ID: ra000005 and hp200271). 
The developers of FE-Project are grateful to Team SCALE for maintaining the SCALE library 
and developers of GFD-Dennou Club providing visualization tools. 
We also thanks Seiya Nishizawa, Hiroaki Miura and Yukio Masumoto 
for their valuable comments and suggestions. 