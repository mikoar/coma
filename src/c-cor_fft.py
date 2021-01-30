#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np

# do kolejnych elementow sekwencji aplikuje funkcje haszujaca o wartosciach +/-1


def hashseq(seq):
    def hashn(w):
        return (hash(w) % 2) * 2 - 1
    return map(hashn, seq)

# liczy korelacje krzyzowa za pomoca fft


def crossfft(a, b):
    N_cross = len(a)+len(b)-1
    N_fft = 2**int(np.log2(N_cross)+1)
    ffta = np.fft.rfft(a, N_fft)
    fftb = np.fft.rfft(b[slice(None, None, -1)], N_fft)
    fftc = ffta*fftb
    return np.round(np.fft.irfft(fftc)[:N_cross]).astype(int)

# liczy znormalizowane wartosci korelacji krzyzowej
# model zaklada, ze argumenty zawieraja tylko wartosci +/-1
# czyli jest nieadekwatny do naszego problemu


def normalized_corr(refer_h, query_h):
    len_refer = len(refer_h)
    len_query = len(query_h)
    corr = crossfft(refer_h, query_h)
    len_corr = len(corr)
    len_min = min(len_refer, len_query)
    len_max = max(len_refer, len_query)

    def n_par(i):
        if i < len_min:
            return i+1+2  # pseudocount
        elif i >= len_max:
            return len_corr-i+2  # pseudocount
        else:
            return len_min+2  # pseudocount
    n_corr = np.array([n_par(i) for i in range(len_corr)])

    def aggreg(hs, len_corr=len_corr):
        ret = []
        s = 1  # pseudocount
        n_hs = len(hs)
        n_ot = len_corr+1 - n_hs
        for i in range(len_corr):
            if i < n_hs and hs[i] == 1:
                s += 1
            if i >= n_ot and hs[i-n_ot] == 1:
                s -= 1
            ret.append(s)
        return np.array(ret)
    f_refer = aggreg(refer_h)
    f_query = aggreg(query_h)
    f_corr = (f_refer*f_query + (n_corr-f_refer) *
              (n_corr-f_query))/(n_corr.astype(float))
    m_corr = 2*f_corr - n_corr
    s_corr = 2*np.sqrt(f_corr*(n_corr-f_corr)/n_corr)
    return (corr - m_corr)/s_corr

# znajduje przesuniecia, dla ktorych znormalizowana wartosc korelacji krzyzowej przekracza thr_corr


def find_shifts(aggr_corrs, thr_corr, len_query):
    shifts = []
    for i, v in enumerate(aggr_corrs):
        if v > thr_corr:
            shifts.append((i-len_query+1, v))
    return shifts

# dla zadanych przesuniec znajduje podobne fragmenty w sekwencjach reference i query o dlugosci > thr_len


def find_matches(shifts, thr_len, reference, query):
    len_query = len(query)
    len_refer = len(reference)
    matches = []
    for s, v in shifts:
        st = max(0, s)
        en = min(len_refer, len_query+s)
        curr_cum = 0
        curr_st = st
        best_cum = 0
        for i in range(st, en):
            if reference[i] == query[i-s]:
                curr_cum += 1
                if curr_cum > best_cum:
                    best_cum = curr_cum
                    best_st = curr_st
                    best_en = i
            else:
                curr_cum -= 1
                if curr_cum < 0:
                    curr_cum = 0
                    curr_st = i+1
                    if best_cum > thr_len:
                        matches.append((best_cum, best_st, best_en+1, s))
                        best_cum = 0
        if best_cum > thr_len:
            matches.append((best_cum, best_st, best_en+1, s))
    return matches

# czysci liste znalezionych dopasowan z zawierajacych nakladajace sie fragmenty query
# usuwane sa dopasowania o nalozeniach >thr_overlap*dlugosc


def filter_overlaps(matches, thr_overlap):
    def nearly_disjoint(m, mf):
        ovlp_st = max(m[1]-m[3], mf[1]-mf[3])
        ovlp_en = min(m[2]-m[3], mf[2]-mf[3])
        return (ovlp_en-ovlp_st) < thr_overlap*(m[2]-m[1])
    filtered = []
    for m in sorted(matches, reverse=True):
        if all(nearly_disjoint(m, mf) for mf in filtered):
            filtered.append(m)
    return filtered


# Przyklad

# parametry
thr_corr = 0.9
thr_len = 5
thr_overlap = 0.5

# sekwencje
reference = 'ala_ma_psa_a_pies_ma_kota'
query = 'ala_ma_kota'

# szukamy podobnych fragmentow
refer_h = hashseq(reference)
query_h = hashseq(query)
norm_corrs = normalized_corr(refer_h, query_h)
shifts = find_shifts(norm_corrs, thr_corr, len(query))
matches = find_matches(shifts, thr_len, reference, query)

# wyniki
print 'Found matches:'
for m in matches:
    print m
    print reference[m[1]: m[2]]
    print query[m[1]-m[3]: m[2]-m[3]]

# po odfiltrowaniu nalozen
filtered = filter_overlaps(matches, 0.5)
print 'Filtered matches:'
for m in filtered:
    print m
    print reference[m[1]: m[2]]
    print query[m[1]-m[3]: m[2]-m[3]]


# %%
