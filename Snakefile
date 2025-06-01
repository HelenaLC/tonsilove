# cosmx-ton

# preamble =========================================

# environment
import re
import pandas as pd
configfile: "config.yaml"
R = config["R"]

# slide metadata, including 
# FOV-to-section mapping
md = "data/raw/metadata.csv"
MD = pd.read_csv(md)
DID = list(set(MD["did"]))
SID = MD["sid"]

SUB = ["bcs", "tcs", "mye", "str"]

# targets ==========================================

# outputs (outs/*.rds)
sce = "data/ref/sce.rds"
mtx = "data/ref/mtx.rds"
mty = "data/ref/mty,{sub}.rds"
raw = "outs/raw-{sid}"
fil = "outs/fil-{sid}.rds"
ccc = "outs/ccc-{sid}.rds"
ist = "outs/ist-{sid}.rds"
lv1 = "outs/lv1-{sid}.rds"
ctx = "outs/ctx-{sid}.rds"
sub = "outs/sub-{sid},{sub}.rds"
jst = "outs/jst-{sid},{sub}.rds"
lv2 = "outs/lv2-{sid},{sub}.rds"
pol = "outs/pol-{sid}.parquet"

pbs = "outs/pbs.rds"
qbs = "outs/qbs-{sub}.rds"

# visuals (plts/*.pdf)
plt = []

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
foo = "code/10-plt__{x}__{y}-{{a}},{{b}},{{p}}.R"
for xy in [
    ["sid", "sid"],
    ["sid", "one_sid_all_sub"],
    ["one_sid_all_sub", "one_sid_all_sub"]]:
    bar = glob_wildcards(expand(foo, x=xy[0], y=xy[1])[0])
    for a,b,p in zip(bar.a, bar.b, bar.p):
        plt += expand(pdf, out1=a, out2=b, plt=p, sid=SID)

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

rule all:
	input:
		pbs,
		sce, mtx,
		#sce, mtx, expand(mty, sub=SUB), 
		expand([raw, fil, ccc], sid=SID),
		plt
		#expand([ist, lv1, ctx], sid=SID),	
		#expand([pbs, qbs], sub=SUB),
		#expand([sub, jst], sid=SID, sub=SUB)

def pool(x): return(x if type(x) == str else ";".join(x))

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
	threads: 10
	priority: 98
	input: 	"code/01-mtx.R", sce
	output:	mtx
	log:    "logs/mtx.Rout"
	shell: '''{R} CMD BATCH\\
	--no-restore --no-save "--args ths={threads}\
	{input[1]} {output}" {input[0]} {log}'''

rule mty:
	threads: 10
	priority: 97
	input: 	"code/02-mty,{sub}.R", sce, mtx
	output:	mty
	log:    "logs/mty,{sub}.Rout"
	shell: '''{R} CMD BATCH\\
	--no-restore --no-save "--args ths={threads}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''

# analysis =========================================

# coercion
for did in DID:
	rule:
		priority: 99
		name: "raw-%s" % did
		input:	"code/00-raw.R", "data/raw/"+did, md
		output:	directory(expand(raw, sid=MD["sid"][MD["did"] == did]))
		log:    expand("logs/raw-{did}.Rout", did=did)
		shell: '''R CMD BATCH\
		--no-restore --no-save "--args wcs={wildcards}\
		{input[1]} {input[2]} {output}" {input[0]} {log}'''

# polygons
def foo(wildcards): return(expand(
	"data/raw/{did}/polygons.csv.gz", 
	did=MD["did"][MD["sid"] == wildcards.sid]))
rule pol:
	priority: 98
	input:	"code/01-pol.R", foo, raw
	output:	pol
	log:    "logs/pol-{sid}.Rout"
	shell: '''R CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''

# filtering
rule fil:
	priority: 98
	input: 	"code/01-fil.R", raw
	output:	fil
	log:    "logs/fil-{sid}.Rout"
	shell: '''{R} CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {output}" {input[0]} {log}'''

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
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''	

# contexts
rule ctx:
    priority: 95
	input:	"code/03-ctx.R",
			a = expand(fil, sid=SID),
			b = expand(ist, sid=SID)
	params:	lambda wc, input: pool(input.a),
			lambda wc, input: pool(input.b)
	output:	expand(ctx, sid=SID)
	log:    "logs/ctx.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{params[0]} {params[1]} {output}" {input[0]} {log}'''	

# subanalysis ======================================

# subsetting
rule sub:
    priority: 95
    input:	"code/04-sub.R", fil, lv1
    output:	expand("outs/sub-{{sid}},{sub}.rds", sub=SUB)
    log:    "logs/sub-{sid}.Rout"
    shell: '''R CMD BATCH\\
    --no-restore --no-save "--args\
    {input[1]} {input[2]} {output}" {input[0]} {log}'''	

# subclustering
rule jst:
	priority: 94
	input:	"code/05-jst.R", mty, 
			a = expand("outs/sub-{sid},{{sub}}.rds", sid=SID)
	params:	lambda wc, input: pool(input.a)
	output:	expand("outs/jst-{sid},{{sub}}.rds", sid=SID)
	log:    "logs/jst-{sub}.Rout"
	shell: '''{R} CMD BATCH --no-restore --no-save "--args\
	{input[1]} {params} {output}" {input[0]} {log}'''

# labeling
rule lv2:
	priority: 93
	input:	"code/06-lv2.R", jst, 
			"meta/labs/lv2,{sub}.json"
	output:	lv2
	log:    "logs/lv2-{sid},{sub}.Rout"
	shell: '''{R} CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''

# profiles	
foo = {
	"one": [
		expand("outs/fil-{sid}.rds", sid=SID),
		expand("outs/lv1-{sid}.rds", sid=SID)],
	"two": [
		expand("outs/sub-{sid},{{sub}}.rds", sid=SID),
		expand("outs/jst-{sid},{{sub}}.rds", sid=SID)]
	}
for bar in foo.keys():
	rule:
		priority: 90
		threads: 10
		name: "pbs-%s" % bar
		input:	"code/00-pbs.R", a=foo[bar][0], b=foo[bar][1]
		params:	lambda wc, input: pool(input.a),
				lambda wc, input: pool(input.b)
		output:	{"one": pbs, "two": qbs}[bar]
		log:    "logs/%s.Rout" % {"one": "pbs", "two": "qbs-{sub}"}[bar]
		shell: '''R CMD BATCH\\
		--no-restore --no-save "--args wcs={wildcards}\
		{params} {output} ths={threads}" {input[0]} {log}'''

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

rule plt_one:
	priority: 99
	input:	"10-plt__one-{out},{plt}.R", "outs/{out}.rds"
	log:	"logs/plt__one-{out},{plt}.Rout"
	output: "plts/{out},{plt}.pdf"
	shell: '''R CMD BATCH\
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
        shell: '''R CMD BATCH\
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
        shell: '''R CMD BATCH\
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
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

# plt by sid
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
        shell: '''R CMD BATCH\
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
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

# def out(by, n=None):
#     o = "out" if n is None else "out"+str(n)
#     if raw:
#     	d = {"all": expand("outs/{{"+o+"}}-{x}/se.rds", x=SID) }
#     else:
# 	    d = {
# 			"one": "outs/{"+o+"}-{x}.rds",
# 			"all": expand("outs/{{"+o+"}}-{x}.rds", x=SID)
# 	    }
#     return(d[by])
#         # "sip": "outs/{"+o+"}-{x}.rds",
#         # "sid": "outs/{"+o+"}-{x}.rds",
#         # "sub": "outs/{"+o+"}-{x}.rds",
#         # "sid_sub": "outs/{"+o+"}-{x}.rds",
#         # "all_sip": expand("outs/{{"+o+"}}-{x}.rds", x=SIP),
#         # "all_sid": expand("outs/{{"+o+"}}-{x}.rds", x=SID),
#         # "all_sub": expand("outs/{{"+o+"}}-{x}.rds", x=SUB),
# 		# "all_raw": expand("outs/{{"+o+"}}-{x}/se.rds", x=SID),
# 		# "all_sid_one_sub": expand("outs/{{"+o+"}}-{y},{{x}}.rds", y=SID),
# 		# "one_sid_all_sub": expand("outs/{{"+o+"}}-{{x}},{y}.rds", y=SUB),
# 		# "all_sid_all_sub": expand("outs/{{"+o+"}}-{x},{y}.rds", x=SID, y=SUB)

# pat = "plt_all-{out},{plt}"
# rule plt_all_one:
# 	priority: 1
# 	input:	"code/10-%s.R" % pat, a = out("all")
# 	params: lambda wc, input: pool(input.a)
# 	log:	"logs/%s.Rout" % pat
# 	output:	"plts/{out},{plt}.pdf"
# 	shell: '''R CMD BATCH\
# 	--no-restore --no-save "--args wcs={wildcards}\
# 	{params} {output}" {input[0]} {log}'''

# pat = "plt_all-{out1},{out2},{plt}"
# rule plt_all_two:
# 	priority: 1
# 	wildcard_constraints: out1 = "^raw"
# 	input:  "code/10-%s.R" % pat,
# 			a = out("all", 1), 
# 			b = out("all", 2)
# 	params: lambda wc, input: pool(input.a),
# 			lambda wc, input: pool(input.b)
# 	log:    "logs/%s.Rout" % pat
# 	output:	"plts/{out1},{out2},{plt}.pdf"
# 	shell: '''R CMD BATCH\
# 	--no-restore --no-save "--args wcs={wildcards}\
# 	{params} {output}" {input[0]} {log}'''

# pat = "plt_all-{out1},{out2},{plt}"
# rule plt_all_raw:
# 	priority: 1
# 	wildcard_constraints: out1 = "raw"
# 	input:  "code/10-%s.R" % pat,
# 			a = out("all", 1, True), 
# 			b = out("all", 2)
# 	params: lambda wc, input: pool(input.a),
# 			lambda wc, input: pool(input.b)
# 	log:    "logs/%s.Rout" % pat
# 	output:	"plts/{out1},{out2},{plt}.pdf"
# 	shell: '''R CMD BATCH\
# 	--no-restore --no-save "--args wcs={wildcards}\
# 	{params} {output}" {input[0]} {log}'''

# pat = "plt_one_one-{out1},{out2},{plt}"
# rule plt_one_one:
# 	input:  "code/10-%s.R" % pat,
# 			a = out("one", 1), 
# 			b = out("one", 2)
# 	params: lambda wc, input: pool(input.a),
# 			lambda wc, input: pool(input.b)
# 	log:    "logs/%s,{x}.Rout" % pat
# 	output: "plts/{out1},{out2},{plt},{x}.pdf"
# 	shell: '''R CMD BATCH\
# 	--no-restore --no-save "--args wcs={wildcards}\
# 	{params} {output}" {input[0]} {log}'''
