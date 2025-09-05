# cosmx-ton

# preamble =========================================

# environment
import re
import json
import pandas as pd
configfile: "config.yaml"
R = config["R"]

# slide metadata, including 
# FOV-to-section mapping
md = "data/raw/metadata.csv"
MD = pd.read_csv(md)
DID = list(set(MD["did"]))
SID = MD["sid"]

SUB = json.load(open("meta/lab/sub.json")).keys()

# regions of interest
ROI = {}
rmv = []
# these exhibit transitions in H&E,
# but not the adjacent section's IF
for sid in SID:
	did = MD["did"][MD["sid"] == sid]
	foo = "imgs/{did}/shapes/{sid}/{{roi}}.pickle"
	foo = expand(foo, did=did, sid=sid)[0]
	foo = glob_wildcards(foo).roi
	foo = [x for x in foo if x not in rmv]
	if sid not in ROI.keys():
		ROI[sid] = {}
	ROI[sid] = foo

# targets ==========================================

# Human Tonsil Atlas-based
# scRNA-seq reference profiles
sce = "data/ref/sce.rds"
mtx = "data/ref/mtx.rds"
mty = "data/ref/mty,{sub}.rds"

# outputs (outs/*.rds)
res = "outs/res-{sid}"
raw = "outs/raw-{sid}"
fil = "outs/fil-{sid}.rds"
pol = "outs/pol-{sid}.parquet"
# regions
roi = "outs/roi-{sid}.rds"
gcs = "outs/gcs-{sid}.rds"
epi = "outs/epi-{sid}.rds"
# downstream
ctx = "outs/ctx-{sid}.rds"
cty = "outs/cty-{sid}.rds"
pro = "outs/pro-{sid}.rds"
ccc = "outs/ccc-{sid}.rds"
# clustering
ist = "outs/ist-{sid}.rds"
lv1 = "outs/lv1-{sid}.rds"
# subclustering
sub = "outs/sub-{sid},{sub}.rds"
jst = "outs/jst-{sid},{sub}.rds"
lv2 = "outs/lv2-{sid},{sub}.rds"
# pseudobulks
pbs = "outs/pbs.rds"
pbt = "outs/pbt.rds"
qbs = "outs/qbs-{sub}.rds"
qbt = "outs/qbt-{sub}.rds"

# visuals (plts/*.pdf)
plt = []
qlt = []

# one
pdf = "plts/{out},{plt}.pdf"
foo = "code/10-plt__one-{a},{p}.R"
bar = glob_wildcards(foo)
for a,p in zip(bar.a, bar.p):
    plt += expand(pdf, out=a, plt=p)

# one by nan
pdf = "plts/{out},{plt},{ids}.pdf"
foo = "code/10-plt__{x}-{{a}},{{p}}.R"
for x in ["sid", "sub"]:
    i = {"sid": SID, "sub": SUB}[x]
    bar = glob_wildcards(expand(foo, x=x)[0])
    for a,p in zip(bar.a, bar.p):
        plt += expand(pdf, out=a, plt=p, ids=i)

# all by nan
pdf = "plts/{out},{plt}.pdf"
foo = "code/10-plt__{x}-{{a}},{{p}}.R"
for x in ["all_sid", "all_sub", "all_sid_all_sub"]:
    bar = glob_wildcards(expand(foo, x=x)[0])
    for a,p in zip(bar.a, bar.p):
        plt += expand(pdf, out=a, plt=p)

# two all by nan
pdf = "plts/{out1},{out2},{plt}.pdf"
foo = "code/10-plt__{x}__{y}-{{a}},{{b}},{{p}}.R"
for xy in [
    ["all_raw", "all_sid"],
	["all_sid", "all_sid"],
	["all_sid", "all_sid_all_sub"]]:
    bar = glob_wildcards(expand(foo, x=xy[0], y=xy[1])[0])
    for a,b,p in zip(bar.a, bar.b, bar.p):
        plt += expand(pdf, out1=a, out2=b, plt=p)

# two by sid
pdf = "plts/{out1},{out2},{plt},{sid}.pdf"
qdf = "qlts/{out1},{out2},{plt},{sid}.pdf"
foo = "code/10-plt__{x}__{y}-{{a}},{{b}},{{p}}.R"
fuu = "code/10-qlt__{x}__{y}-{{a}},{{b}},{{p}}.R"
for xy in [
	["sid", "sid"],
	["sid", "one_sid_all_sub"]]:
	bar = glob_wildcards(expand(foo, x=xy[0], y=xy[1])[0])
	for a,b,p in zip(bar.a, bar.b, bar.p):
		plt += expand(pdf, out1=a, out2=b, plt=p, sid=SID)
	bar = glob_wildcards(expand(fuu, x=xy[0], y=xy[1])[0])
	for a,b,p in zip(bar.a, bar.b, bar.p):
		qlt += expand(qdf, out1=a, out2=b, plt=p, sid=SID)

# one by sub
pdf = "plts/{out1},{plt},{{sub}}.pdf"
foo = glob_wildcards("code/10-plt__all_sid_one_sub-{x},{y}.R")
for x,y in zip(foo.x, foo.y):
	plt += expand(expand(pdf, out1=x, plt=y), sub=SUB)

# # one by sid-sub
# pdf = "plts/{out1},{plt},{{sid}},{{sub}}.pdf"
# foo = glob_wildcards("code/10-plt__sid_sub-{x},{y}.R")
# for x,y in zip(foo.x, foo.y):
# 	plt += expand(expand(pdf, out1=x, plt=y), sid=SID, sub="epi")

# two by sid-sub
pdf = "plts/{out1},{out2},{plt},{{sid}},{{sub}}.pdf"
foo = glob_wildcards("code/10-plt__sid_sub__sid_sub-{x},{y},{z}.R")
for x,y,z in zip(foo.x, foo.y, foo.z):
	plt += expand(expand(pdf, out1=x, out2=y, plt=z), sid=SID, sub=SUB)

# for ex in ["ccc"]:
# 	pat = re.compile(r"^((?!%s).)*$" % ex)
# 	plt = [p for p in plt if pat.match(p)]

rule all:
	input:
	    # visuals
		plt, qlt,
		# reference
		sce, mtx, expand(mty, sub=SUB), 
		# setup
		expand([raw, fil, pol], sid=SID),
		# clustering
		expand([ist, lv1], sid=SID),
		# subclustering
		expand([jst, lv2], sid=SID, sub=SUB),
		# pseudobulks
		pbs, pbt, expand([qbs, qbt], sub=SUB),
		# regions
		expand([roi, epi, gcs], sid=SID),
		# downstream
		expand([ctx, cty, ccc], sid=SID),
		# collection
		expand([res], sid=SID)

def pool(x): return(x if type(x) == str else ";".join(x))

# write session info to .txt file
rule inf:
    input:	"code/_inf.R"
    output:	"sess_info.txt"
    log:	"logs/inf.Rout" 
    shell:	'''R CMD BATCH\\
    --no-restore --no-saves "--args\\
    {output}" {input} {log}'''

# reference ========================================

# loading
rule sce:
	priority: 99
	input: 	"code/00-sce.R",
			expand(raw, sid=SID[0])
	output:	sce
	log:    "logs/sce.Rout"
	shell: '''{R} CMD BATCH\\
	--no-restore --no-save "--args\
	{input[1]} {output}" {input[0]} {log}'''

# aggregation
rule mtx:
	threads: 5
	priority: 98
	input: 	"code/01-mtx.R", sce
	output:	mtx
	log:    "logs/mtx.Rout"
	shell: '''{R} CMD BATCH\\
	--no-restore --no-save "--args ths={threads}\
	{input[1]} {output}" {input[0]} {log}'''

rule mty:
	threads: 5
	priority: 97
	input: 	"code/02-mty,{sub}.R", pbt, sce,
			"meta/lab/sub.json"
	output:	mty
	log:    "logs/mty,{sub}.Rout"
	shell: '''{R} CMD BATCH\\
	--no-restore --no-save "--args ths={threads} wcs={wildcards}\
	{input[1]} {input[2]} {input[3]} {output}" {input[0]} {log}'''

# setup ============================================

# coercion
for did in DID:
	rule:
		priority: 99
		name: "raw-%s" % did
		input:	"code/00-raw.R", "data/raw/"+did, md
		output:	directory(expand(raw, sid=MD["sid"][MD["did"] == did]))
		log:    expand("logs/raw-{did}.Rout", did=did)
		shell: '''{R} CMD BATCH\
		--no-restore --no-save "--args wcs={wildcards}\
		{input[1]} {input[2]} {output}" {input[0]} {log}'''

# regions
def foo(wildcards): return(expand(
	"imgs/{did}/shapes/{sid}/{roi}",
	did=MD["did"][MD["sid"] == wildcards.sid],
	sid=wildcards.sid, roi=ROI[wildcards.sid]))
rule roi:
	input:	"code/03-roi.R", fil, x = foo
	params: lambda wc, input: pool(input.x)
	output:	roi
	log:	"logs/roi-{sid}.Rout"
	shell: '''{R} CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {params} {output[0]}" {input[0]} {log}'''

# analysis =========================================

# filtering
rule fil:
	priority: 98
	input: 	"code/01-fil.R", raw
	output:	fil
	log:    "logs/fil-{sid}.Rout"
	shell: '''{R} CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {output}" {input[0]} {log}'''

# polygons
def pos_csv(wildcards): return(expand(
	"data/raw/{did}/positions.csv.gz", 
	did=MD["did"][MD["sid"] == wildcards.sid]))
def pol_csv(wildcards): return(expand(
	"data/raw/{did}/polygons.csv.gz", 
	did=MD["did"][MD["sid"] == wildcards.sid]))
rule pol:
	priority: 98
	input:	"code/02-pol.R", fil, pos_csv, pol_csv
	output:	pol
	log:    "logs/pol-{sid}.Rout"
	shell: '''{R} CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {input[3]} {output}" {input[0]} {log}'''

# processing
rule pro:
	threads: 10
	priority: 97
	input: 	"code/02-pro.R", fil, mtx
	output:	pro
	log:    "logs/fil-{sid}.Rout"
	shell: '''{R} CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''

# communication
rule ccc:
	threads: 10
	priority: 97
	input: 	"code/02-ccc.R", fil
	output:	ccc
	log:    "logs/ccc-{sid}.Rout"
	shell: '''{R} CMD BATCH\
	--no-restore --no-save "--args ths={threads}\
	{input[1]} {output}" {input[0]} {log}'''

# clustering
rule ist:
	priority: 97
	input: 	"code/02-ist.R", mtx, a = expand(fil, sid=SID)
	params: lambda wc, input: pool(input.a)
	output:	expand(ist, sid=SID)
	log:    "logs/ist.Rout"
	shell: '''{R} CMD BATCH --no-restore --no-save "--args\
	{input[1]} {params} {output}" {input[0]} {log}'''

# labelling
rule lv1:
	priority: 96
	input:	"code/00-lab.R", ist,
			"meta/lab/lv1.json"
	output:	lv1
	log:    "logs/lv1-{sid}.Rout"
	shell: '''{R} CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''	

# using 'sosta', programmatically identify 
# spatial features corresponding to GCs &
# compute cell-wise distance to mantle midpoint
rule gcs:
	priority: 98
	input: 	"code/03-gcs.R", fil, lv1
	output:	gcs
	log:    "logs/gcs-{sid}.Rout"
	shell: '''{R} CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''

# get epithelial sturctures
rule epi:
	priority: 98
	input: 	"code/03-epi.R", fil, lv1
	output:	epi
	log:    "logs/epi-{sid}.Rout"
	shell: '''{R} CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''

# contexts
rule ctx:
	priority: 95
	input:	"code/03-ctx.R",
			a = expand(fil, sid=SID),
			b = expand(lv1, sid=SID)
	params:	lambda wc, input: pool(input.a),
			lambda wc, input: pool(input.b)
	output:	expand(ctx, sid=SID)
	log:    "logs/ctx.Rout"
	shell: '''{R} CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{params[0]} {params[1]} {output}" {input[0]} {log}'''

# merging
rule cty:
	priority: 94
	input:	"code/03-cty.R",
			a = expand(ctx, sid=SID)
	params:	lambda wc, input: pool(input.a)
	output:	expand(cty, sid=SID)
	log:    "logs/cty.Rout"
	shell: '''{R} CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{params[0]} {output}" {input[0]} {log}'''

# subanalysis ======================================

# subsetting
rule sub:
    priority: 95
    input:	"code/04-sub.R", "meta/lab/sub.json", fil, lv1
    output:	expand("outs/sub-{{sid}},{sub}.rds", sub=SUB)
    log:    "logs/sub-{sid}.Rout"
    shell: '''{R} CMD BATCH\\
    --no-restore --no-save "--args\
    {input[1]} {input[2]} {input[3]}\
    {output}" {input[0]} {log}'''	

# subclustering
rule jst:
	priority: 94
	input:	"code/05-jst.R", mty, 
			a = expand("outs/sub-{sid},{{sub}}.rds", sid=SID)
	params:	lambda wc, input: pool(input.a)
	output:	expand("outs/jst-{sid},{{sub}}.rds", sid=SID)
	log:    "logs/jst-{sub}.Rout"
	shell: '''{R} CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {params} {output}" {input[0]} {log}'''

# labeling
rule lv2:
	priority: 93
	input:	"code/00-lab.R", jst, 
			"meta/lab/lv2,{sub}.json"
	output:	lv2
	log:    "logs/lv2-{sid},{sub}.Rout"
	shell: '''{R} CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''

# profiles	
la = [
	"outs/fil-{sid}.rds",
	"outs/sub-{sid},{{sub}}.rds"]
le = {
	"pbs": "outs/ist-{sid}.rds",
	"pbt": "outs/lv1-{sid}.rds",
	"qbs": "outs/jst-{sid},{{sub}}.rds",
	"qbt": "outs/lv2-{sid},{{sub}}.rds"}
lu = {
	"pbs": "logs/pbs.Rout", 
	"pbt": "logs/pbt.Rout", 
	"qbs": "logs/qbs-{sub}.Rout", 
	"qbt": "logs/qbt-{sub}.Rout"}
for x in ["pbs", "pbt", "qbs", "qbt"]:
	a = la[0] if x[0] == "p" else la[1]
	rule:
		priority: 99
		threads: 10
		name: x
		input:	"code/00-pbs.R", 
				a=expand(a, sid=SID), 
				b=expand(le[x], sid=SID)
		params:	lambda wc, input: pool(input.a),
				lambda wc, input: pool(input.b)
		output:	globals()[x]
		log:    lu[x]
		shell: '''{R} CMD BATCH\\
		--no-restore --no-save "--args wcs={wildcards}\
		{params} {output} ths={threads}" {input[0]} {log}'''

# collection
rule res:
    priority: 10
    input:  "code/09-res.R", "outs/raw-{sid}",
            expand("outs/{out}-{{sid}},{sub}.rds", out=["jst", "lv2"], sub=SUB),
            expand("outs/{out}-{{sid}}.rds", out=["fil", "ist", "lv1", "epi", "gcs", "cty", "ccc"])
    output: directory("outs/res-{sid}")
    log:    "logs/res-{sid}.Rout"
    shell: '''R CMD BATCH\\
    --no-restore --no-save "--args wcs={wildcards}\
    {input} {output}" {input[0]} {log}'''

# visualization ====================================

def out(by, n=None):
    o = "out" if n is None else "out"+str(n)
    d = {
        "sid": "outs/{"+o+"}-{x}.rds",
        "sub": "outs/{"+o+"}-{x},{y}.rds",
        "sid_sub": "outs/{"+o+"}-{x},{y}.rds",
        "all_sid": expand("outs/{{"+o+"}}-{x}.rds", x=SID),
        "all_sub": expand("outs/{{"+o+"}}-{x}.rds", x=SUB),
        "all_raw": expand("outs/{{"+o+"}}-{x}/se.rds", x=SID),
        "all_sid_one_sub": expand("outs/{{"+o+"}}-{y},{{x}}.rds", y=SID),
        "one_sid_all_sub": expand("outs/{{"+o+"}}-{{x}},{y}.rds", y=SUB),
        "all_sid_all_sub": expand("outs/{{"+o+"}}-{x},{y}.rds", x=SID, y=SUB)
    }
    return(d[by])

rule plt__one:
	priority: 99
	input:	"10-plt__one-{out},{plt}.R", "outs/{out}.rds"
	log:	"logs/plt__one-{out},{plt}.Rout"
	output: "plts/{out},{plt}.pdf"
	shell: '''{R} CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {output}" {input[0]} {log}'''

pat = "plt__{by}-{{out}},{{plt}}"
for by in ["sid", "all_sid_one_sub"]:
    rule:
        name:   "plt__%s" % by
        priority: 99
        input:  expand("code/10-"+pat+".R", by=by), a = out(by)
        params: lambda wc, input: pool(input.a)
        log:    expand("logs/"+pat+",{{x}}.Rout", by=by) 
        output: "plts/{out},{plt},{x}.pdf"
        shell: '''{R} CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

pat = "plt__{by}-{{out}},{{plt}}"
for by in ["all_sid", "all_sub", "all_sid_all_sub"]:
    rule:
        name:   "plt__%s" % by
        priority: 99
        input:	expand("code/10-"+pat+".R", by=by), a = out(by)
        params: lambda wc, input: pool(input.a)
        log:	expand("logs/"+pat+".Rout", by=by)
        output:	"plts/{out},{plt}.pdf"
        shell: '''{R} CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

pat = "plt__{by1}__{by2}-{{out1}},{{out2}},{{plt}}"
for by in [
    ["all_raw","all_sid"],
	["all_sid","all_sid"],
	["all_sid","all_sid_all_sub"]]:
    rule:
        name:   "plt__%s__%s" % (by[0], by[1])
        priority: 99
        input:  expand("code/10-"+pat+".R", by1=by[0], by2=by[1]),
                a = out(by[0], 1), b = out(by[1], 2)
        params: lambda wc, input: pool(input.a),
                lambda wc, input: pool(input.b)
        log:    expand("logs/"+pat+".Rout", by1=by[0], by2=by[1])
        output:	"plts/{out1},{out2},{plt}.pdf"
        shell: '''{R} CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

pat = "plt__{by1}__{by2}-{{out1}},{{out2}},{{plt}}"
for by in [
    ["sid","sid"],
    ["sid","one_sid_all_sub"],
    ["one_sid_all_sub","one_sid_all_sub"]]:
    rule:
        name:   "plt__%s__%s" % (by[0], by[1])
        priority: 99
        input:  expand("code/10-"+pat+".R", by1=by[0], by2=by[1]),
                a = out(by[0], 1), b = out(by[1], 2)
        params: lambda wc, input: pool(input.a),
                lambda wc, input: pool(input.b)
        log:    expand("logs/%s,{{x}}.Rout" % pat, by1=by[0], by2=by[1])
        output: "plts/{out1},{out2},{plt},{x}.pdf"
        shell: '''{R} CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

pat = "plt__{by1}__{by2}-{{out1}},{{out2}},{{plt}}"
for by in [
	["sid_sub","sid_sub"]]:
    rule:
        name:   "plt__%s__%s" % (by[0], by[1])
        priority: 99
        input:  expand("code/10-%s.R" % pat, by1=by[0], by2=by[1]),
                a = out(by[0], 1), b = out(by[1], 2)
        params: lambda wc, input: pool(input.a),
                lambda wc, input: pool(input.b)
        log:    expand("logs/%s,{{x}},{{y}}.Rout" % pat, by1=by[0], by2=by[1])
        output: "plts/{out1},{out2},{plt},{x},{y}.pdf"
        shell: '''{R} CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

pat = "qlt__{by1}__{by2}-{{out1}},{{out2}},{{plt}}"
for by in [
    ["sid","sid"],
    ["sid","one_sid_all_sub"]]:
    rule:
        name:   "qlt__%s__%s" % (by[0], by[1])
        priority: 99
        input:  expand("code/10-"+pat+".R", by1=by[0], by2=by[1]),
                "outs/pol-{x}.parquet", 
                a = out(by[0], 1),
                b = out(by[1], 2),
        params: lambda wc, input: pool(input.a),
                lambda wc, input: pool(input.b)
        log:    expand("logs/%s,{{x}}.Rout" % pat, by1=by[0], by2=by[1])
        output: "qlts/{out1},{out2},{plt},{x}.pdf"
        shell: '''{R} CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {input[1]} {output}" {input[0]} {log}'''
