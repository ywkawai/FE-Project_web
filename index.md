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
This project provides a library and sample programs for the discontinuous Galerkin methods (DGM). Futheremore, we provide an atmospheric model with a regional/global dynamical core based on DGM, SCALE-DG. 
</p>

![Top Image]({{ '/assets/images/FE-project_top.png' | relative_url }})


# Contributors

- Main developers
  - [Yuta Kawai](https://researchmap.jp/ykawai1988/?lang=english)<sup>*1</sup> 
- Project design and Scientific decision 
  - Yuta Kawai<sup>*1</sup>  and Hirofumi Tomita<sup>*1</sup>  
- Investigation and improvement of computational performance 
  - Xuanzhengbo Ren<sup>*2</sup>

Affiliation: *1: [RIKEN Center for Computational Science](http://www.r-ccs.riken.jp/en/), Kobe in Japan; *2: Nagoya University, Nagoya in Japan

# Acknowledgements

This project is supported by 
[the Transformaive Research Areas B: DNA Climate Science](https://dna-climate.org/) (MEXT KAKENHI Grant Number JP20H05731), 
[Moonshot Goal8 Realization of a society safe from the threat of extreme winds and rains by controlling and modifying the weather by 2050](https://www.jst.go.jp/moonshot/program/goal8/) ([Development of an atmospheric simulation model for probability estimation for local atmospheric phenomena](https://moonshot8-modeldev.riken.jp)), and JST AIP Grant Number JPMJCR19U2. 
The model development and numerical experiments are
performed using supercomputers (Oackbridge-CX and Wisteria) at the University of Tokyo and Fugaku at RIKEN (Project ID: ra000005, hp200271, hp230278). 
The developers of FE-Project are grateful to Team SCALE for maintaining the SCALE library 
and developers of GFD-Dennou Club providing visualization tools. 
We also thank Dr. Seiya Nishizawa, Hiroaki Miura, Keiichi Ishioka, and Yukio Masumoto 
for their valuable comments and suggestions. 