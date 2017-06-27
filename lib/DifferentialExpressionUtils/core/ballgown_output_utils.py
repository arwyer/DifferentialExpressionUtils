# This downloads and unzips in a separate directory the stringtie output from each expression object
# in the given expression set, and verifies that it has the correct files for ballgown to run
# (*.ctab).  Each directory name is prefixed by the given dir_prefix so that the ballgown call can
# use it to specify the input directories to ballgown
#
# This also returns information extracted from the ExpressionSet object which is useful for further
# workspace storage (genome_id, etc )
#
def get_info_and_download_for_ballgown( logger, ws_client, hs, ws_id, urls, directory, dir_prefix, expressionset_id, token ):
        try:
           #expression_set = ws_client.get_objects( [ { 'name'      : expressionset_id,
           #                                            'workspace' : ws_id }
           #                                        ] )[0]
           expression_set = script_util.ws_get_obj(logger, ws_client, ws_id, expressionset_id)[0]

        except Exception,e:
           raise Exception( "".join(traceback.format_exc() ))
        ### Getting all the set ids and genome_id
        output_obj = {}
        #expression_set_info = ws_client.get_object_info_new({"objects": [{'name' : expressionset_id, 'workspace': ws_id}]})[0]
        #output_obj['expressionset_id'] =  str(expression_set_info[6]) + '/' + str(expression_set_info[0]) + '/' + str(expression_set_info[4])
        output_obj['expressionset_id'] = script_util.ws_get_ref(logger, ws_client, ws_id, expressionset_id)
        #output_obj['sampleset_id'] =  expression_set['data']['sampleset_id']
        #output_obj['alignmentset_id'] = expression_set['data']['alignmentSet_id']
        #output_obj['genome_id'] = expression_set['data']['genome_id']
        #output_obj['genome_name'] = ws_client.get_object_info([{"ref" :output_obj['genome_id']}],includeMetadata=None)[0][1]
        #ws_gtf = output_obj['genome_name']+"_GTF_Annotation"

        # QUESTION:  whats the reason for getting this from object_info and not the object
        output_obj['genome_id']             = expression_set['data']['genome_id']
        output_obj['alignmentSet_id']       = expression_set['data']['alignmentSet_id']
        output_obj['sampleset_id']          = expression_set['data']['sampleset_id']
        output_obj['sample_expression_ids'] = expression_set['data']['sample_expression_ids']
        ### Check if GTF object exists in the workspace pull the gtf
        #ws_gtf = output_obj['genome_name']+"_GTF"
        #gtf_file = script_util.check_and_download_existing_handle_obj(logger,ws_client,urls,ws_id,ws_gtf,"KBaseRNASeq.GFFAnnotation",directory,token)
        #print 'GTF file is  ' +  gtf_file
        #if gtf_file is None:
        #     create_gtf_annotation_from_genome(logger,ws_client,hs,urls,ws_id,output_obj['genome_id'],output_obj['genome_name'],directory,token)
        #output_obj['gtf_file'] = gtf_file
        ### Getting the expression objects and alignment objects
        m_expr_ids = expression_set['data']['mapped_expression_ids']
        if len(m_expr_ids)  < 2:
           raise ValueError("Error the ExpressionSet object has less than 2 expression samples. Kindly check your reads files and repeat the previous step (Cufflinks)")
        output_obj['labels'] = []   # not sure why this is still needed - remove if not
        output_obj['conditions'] = []
        #output_obj['alignments'] = []
        output_obj['subdirs'] = []
        counter = 0
        #assembly_file = os.path.join( directory, "assembly_gtf.txt" )
        #list_file = open( assembly_file, 'w' )
        for i in m_expr_ids:
            for a_id, e_id in i.items():
                #files = {}
                #a_obj,e_obj = ws_client.get_objects(
                #                     [{'ref' : a_id},{'ref': e_id}])
                #e_obj= ws_client.get_objects( [ {'ref' : e_id} ] )[0]
                e_obj = script_util.ws_get_obj( logger, ws_client, ws_id, e_id )[0]

                zipfile = e_obj['data']['file']['file_name']
                shock_id = e_obj['data']['file']['id']
                ### Get the condition name, replicate_id , shock_id and shock_filename
                condition = e_obj['data']['condition']                                            # Question: can we use this for group!?
                if 'replicate_id' in e_obj['data'] : replicate_id = e_obj['data']['replicate_id']
                #files[a_obj['data']['file']['file_name']] = a_obj['data']['file']['id']
                #files[e_obj['data']['file']['file_name']] = e_obj['data']['file']['id']
                if not condition in output_obj['labels']:
                    output_obj['labels'].append( condition )
                else:
                    counter += 1 #### comment it when replicate_id is available from methods
                output_obj['conditions'].append( condition )

                #subdir = os.path.join( directory, dir_prefix + "_" + condition + "_" + str(counter) ) ### Comment this line when replicate_id is available from the methods
                # Fix this - we're getting expression object name from the zip file
                if ( zipfile[-4:] != '.zip' ):
                    raise  Exception( "zip file {0} doesn't seem to have .zip extention, can't form a subdirectory name confidently" )

                subdir = os.path.join( directory, zipfile[0:-4] )
                logger.info( "subdir is {0}".format( subdir ) )
                output_obj['subdirs'].append( subdir )
                if not os.path.exists( subdir ):
                    os.makedirs( subdir )
                try:
                    script_util.download_file_from_shock( logger = logger,
                                                          shock_service_url = urls['shock_service_url'],
                                                          shock_id = shock_id,
                                                          filename = zipfile,
                                                          directory = subdir,
                                                          token = token )
                    #script_util.download_shock_files( logger, urls['shock_service_url'], s_path, files, token )
                except Exception,e:
                    raise Exception( "Unable to download shock file, {0}".format(e))
                try:
                    script_util.unzip_files( logger,
                                             os.path.join( subdir, zipfile ),
                                             subdir )
                    logger.info( "listing of {0}".format( subdir ))
                    logger.info( os.listdir(subdir) )
                    #script_util.unzip_files(logger,os.path.join(s_path,e_obj['data']['file']['file_name']),s_path)
                    #e_file_path =  os.path.join(s_path,"transcripts.gtf")
                    #a_file_path = os.path.join( s_path, "accepted_hits.bam" )
                    #if os.path.exists(a_file_path) :
                    #         print a_file_path
                    #         output_obj['alignments'].append(a_file_path)
                    #if os.path.exists(e_file_path) : list_file.write("{0}\n".format(e_file_path))
                except Exception, e:
                    raise Exception("".join(traceback.format_exc()))

                for f in [ 'e2t.ctab', 'e_data.ctab', 'i2t.ctab', 'i_data.ctab', 't_data.ctab' ]:
                    fullpath = os.path.join( subdir, f )
                    if ( not os.path.isfile( fullpath ) ):
                        raise Exception( "error: ballgown input file {0} not found. Can't proceed".format( fullpath ) )
                    # check here to see if introns is empty.  Currently (1 Feb 2017) ballgown can't handle this
                    if ( f == 'i2t.ctab' ):
                        if ( number_of_lines_in_file( fullpath ) < 2 ):
                            raise Exception( "error: ballgown does not yet support prokarytic data (intron count appears to be zero)" )

        #list_file.close()
        #output_obj['gtf_list_file'] = assembly_file
        logger.info( output_obj )
        return output_obj


def number_of_lines_in_file(fname):
    n = -1
    with open(fname) as f:
        for i, l in enumerate(f):
            n = n + 1
    return n + 1


def create_sample_dir_group_file(logger,
                                 ws_client,
                                 ws_id,
                                 subdir_list,
                                 condition_list,
                                 sample_dir_group_file
                                 ):
    # QUESTION:  is it possible that duplicate subdirs can appear here?
    #            and with same or differen
    # ballgown requires numeric identifiers for experimental condition groups,
    # so first make a table mapping condition names

    logger.info("new create_sample_dir_group_file:")
    logger.info(pformat(condition_list))
    logger.info(pformat(subdir_list))

    ngroups = 0
    group_name_indices = {}
    group_counts = {}

    for group in condition_list:
        if not group in group_name_indices:
            group_name_indices[group] = ngroups
            ngroups = ngroups + 1
        if not group in group_counts:
            group_counts[group] = 1
        else:
            group_counts[group] = group_counts[group] + 1

    # checks for proper ballgown execution:
    if ngroups < 2:
        raise Exception("At least two condition groups are needed for this analysis")
    for group in condition_list:
        if group_counts[group] < 2:
            raise Exception("condition group {0} has less than 2 members; ballgown will not run".format(group))

    # write the file

    try:
        f = open(sample_dir_group_file, "w")
    except Exception:
        raise Exception("Can't open file {0} for writing {1}".format(sample_dir_group_file, traceback.format_exc()))

    for subdir, group in zip(subdir_list, condition_list):
        f.write("{0}  {1}\n".format(subdir, group_name_indices[group]))
    f.close()

    return

    def run_ballgown_diff_exp(logger,
                              rscripts_dir,
                              directory,
                              sample_dir_group_table_file,
                              ballgown_output_dir,
                              output_csv,
                              volcano_plot_file
                              ):

        # sample_group_table is a listing of output Stringtie subdirectories,
        # (full path specification) paired with group label (0 or 1), ie
        #    /path/WT_rep1_stringtie    0
        #    /path/WT_rep2_stringtie    0
        #    /path/EXP_rep1_stringtie   1
        #    /path/EXP_rep2_stringtie   1
        #  (order doesn't matter, but the directory-group correspondance does)

        #  1) Need to make a generic "Run rscript program"
        #  2) is this the best way to run the R script (subprocess Popen?)


        # Make call to execute the system.

        rcmd_list = ['Rscript', os.path.join(rscripts_dir, 'ballgown_fpkmgenematrix.R'),
                     '--sample_dir_group_table', sample_dir_group_table_file,
                     '--output_dir', ballgown_output_dir,
                     '--output_csvfile', output_csv,
                     '--volcano_plot_file', volcano_plot_file
                     ]
        rcmd_str = " ".join(str(x) for x in rcmd_list)
        logger.info("rcmd_string is {0}".format(rcmd_str))
        openedprocess = subprocess.Popen(rcmd_str, shell=True)  # , stdout=subprocess.PIPE )
        openedprocess.wait()
        # Make sure the openedprocess.returncode is zero (0)
        if openedprocess.returncode != 0:
            logger.info("R script did not return normally, return code - "
                        + str(openedprocess.returncode))
            raise Exception("Rscript failure")


# reads csv diff expr matrix file from Ballgown and returns as a
# dictionary of rows with the gene as key.  Each key gives a row of
# length three corresponding to fold_change, pval and qval in string form
# - can include 'NA'
#
def load_diff_expr_matrix(ballgown_output_dir, output_csv):
    diff_matrix_file = os.path.join(ballgown_output_dir, output_csv)

    if not os.path.isfile(diff_matrix_file):
        raise Exception("differential expression matrix csvfile {0} doesn't exist!".format(diff_matrix_file))

    n = 0
    dm = {}
    with  open(diff_matrix_file, "r") as csv_file:
        csv_rows = csv.reader(csv_file, delimiter="\t", quotechar='"')
        for row in csv_rows:
            n = n + 1
            if (n == 1):
                if (row != ['id', 'fc', 'pval', 'qval']):
                    raise Exception("did not get expected column heading from {0}".format(diff_matrix_file))
            else:
                if (len(row) != 4):
                    raise Exception("did not get 4 elements in row {0} of csv file {1} ".format(n, diff_matrix_file))
                key = row[0]
                # put in checks for NA or numeric for row[1] through 4
                if (key in dm):
                    raise Exception("duplicate key {0} in row {1} of csv file {2} ".format(key, n, diff_matrix_file))
                dm[key] = row[1:5]

    return dm


# this takes the output_csv generated by run_ballgown_diff_exp() and loads it into
# an RNASeqDifferentialExpression named output_object_name

def load_ballgown_output_into_ws(logger,
                                 ws_id,
                                 ws_client,
                                 hs_client,
                                 token,
                                 directory,
                                 ballgown_output_dir,
                                 tool_used,
                                 tool_version,
                                 sample_ids,
                                 conditions,
                                 genome_id,
                                 expressionset_id,
                                 alignmentset_id,
                                 sampleset_id,
                                 output_object_name
                                 ):
    logger.info("Zipping ballgown output")
    zip_file_path = os.path.join(directory, "{0}.zip".format(output_object_name))
    try:
        script_util.zip_files(logger, ballgown_output_dir, zip_file_path)
    except Exception, e:
        raise Exception("Error creating zip file of ballgown output")

    logger.info("Creating handle from ballgown output zip file")
    try:
        handle = hs_client.upload(zip_file_path)
    except Exception, e:
        raise Exception(
            "Error uploading ballgown output zip file to handle service: {0}".format(" ".join(traceback.print_exc())))

    # create object

    de_obj = {"tool_used": tool_used,
              "tool_version": tool_version,
              "sample_ids": sample_ids,
              "condition": conditions,
              "genome_id": genome_id,
              "expressionSet_id": expressionset_id,
              "alignmentSet_id": alignmentset_id,
              "sampleset_id": sampleset_id,
              "file": handle
              }

    objs_save_data = ws_client.save_objects(
        {"workspace": ws_id,
         "objects": [
             {
                 "type": "KBaseRNASeq.RNASeqDifferentialExpression",
                 "data": de_obj,
                 "name": output_object_name
             }
         ]
         }
    )
    return objs_save_data[0]


#  returns a list of gene names from the keys of diff_expr_matrix
#  assumed AND logic between all the tests.   The gene_names are
#  ordered by decreasing fold-change

def filter_genes_diff_expr_matrix(diff_expr_matrix,
                                  scale_type,  # "linear", "log2+1", "log10+1"
                                  qval_cutoff,
                                  fold_change_cutoff,  # says this is log2 but we should make it match scale_type
                                  max_genes):
    # verify that qval_cutoff, fold change cutoff is None or numeric.
    # verify scale_type is an allowed value

    # !!!!! figure out what to do with scale conversion

    if max_genes == None:
        max_genes = sys.maxint

    selected = []
    ngenes = 0

    logbase = 0  # by default this indicates linear
    if (scale_type.lower()[0:4] == "log2"):
        logbase = 2
    elif (scale_type.lower()[0:5] == "log10"):
        logbase = 10

    # iterate through keys (genes) in the diff_expr_matrix, sorted by the first value of each row (fc)
    # (decreasing)

    for gene, v in sorted(diff_expr_matrix.iteritems(), key=lambda (k, v): (convert_NA_low(v[0]), k), reverse=True):

        fc, pval, qval = diff_expr_matrix[gene]

        if logbase > 0 and fc != 'NA' and fc != 'Nan':
            fc = make_numeric(fc, "fc, gene {0}, about to take log{1}".format(gene, logbase))
            try:
                fc = math.log(fc + 1, logbase)
            except:
                raise Exception("unable to take log{0} of fold change value {1}".format(logbase, fc))

        # should NA, NaN fc, qval automatically cause exclusion?

        if (qval_cutoff != None):  # user wants to apply qval_cutoff
            if qval == 'NA' or qval == "Nan":  # bad values automatically excluded
                continue
            q = make_numeric(qval, "qval, gene {0}".format(gene))
            if (q < qval_cutoff):
                continue

        if (fold_change_cutoff != None):  # user wants to apply fold_change_cutoff
            if fc == 'NA' or fc == "Nan":  # bad values automatically excluded
                continue
            f = make_numeric(fc, "fc, gene {0}".format(gene))
            if (f < fold_change_cutoff):
                continue

        ngenes = ngenes + 1
        if ngenes > max_genes:
            break

        # if we got here, it made the cut, so add the gene to the list
        selected.append(gene)

    return selected


# This weeds out all the data (rows, mappings) in given expression matrix object
# for genes not in the given selected_list, and returns a new expression matrix object
# The rows and values presevere the order of the input selected_list

def filter_expr_matrix_object(emo, selected_list):
    fmo = {}
    fmo["type"] = emo["type"]
    fmo["scale"] = emo["scale"]
    fmo["genome_ref"] = emo["genome_ref"]

    fmo["data"] = {}
    fmo["data"]["col_ids"] = emo["data"]["col_ids"]
    fmo["data"]["row_ids"] = []
    fmo["data"]["values"] = []
    fmo["feature_mapping"] = {}

    nrows = len(emo["data"]["row_ids"])
    if nrows != len(emo["data"]["values"]):
        raise Exception("filtering expression matrix:  row count mismatch in expression matrix")

    # make index from gene to row
    gindex = {}
    for i in xrange(nrows):
        gindex[emo["data"]["row_ids"][i]] = i

    for gene in selected_list:
        try:
            i = gindex[gene]
        except:
            raise Exception("gene {0} from differential expression not found in expression matrix")
        fmo["data"]["row_ids"].append(emo["data"]["row_ids"][i])
        fmo["data"]["values"].append(emo["data"]["values"][i])
        fmo["feature_mapping"][gene] = emo["feature_mapping"][gene]

    return fmo


def create_and_save_volcano_plot_report(logger,
                                        ws_client,
                                        ws_id,
                                        callback_url,
                                        token,
                                        ballgown_output_dir,
                                        volcano_plot_file,
                                        de_obj_ref,
                                        em_obj_ref,
                                        report_obj_name):
    logger.info("in create_and_save_volcano_plot_report, callback url is {0}, plot file is {1}".format(callback_url,
                                                                                                       volcano_plot_file))
    # probably DTU needs to be called here?
    volcano_file_path = os.path.join(ballgown_output_dir, volcano_plot_file)
    #  How best to zip this
    # if more than one, use script_util.zip_files
    ##logger.info( "zipping volcano plot file")
    ##image_zip_file = volcano_plot_file + ".zip"   # omit path
    ##image_zip_path = os.path.join( ballgown_output_dir, image_zip_file)
    ##
    ##with ZipFile( image_zip_path, 'w', allowZip64=True) as izip:
    ##    izip.write( volcano_file_path, volcano_plot_file )
    ##
    ##logger.info( "making shock handle for zipped ")
    # image_zip_shock_ret = script_util.upload_file_to_shock( logger, image_zip_path )
    # image_zip_shock_ret = script_util.upload_file_to_shock( logger, image_zip_path )

    if os.path.exists(volcano_file_path):
        volcano_file_shock_ret = script_util.upload_file_to_shock(logger, volcano_file_path)
        logger.info(pformat(volcano_file_shock_ret))
    else:
        logger.info("no volcano file {0} found - skipping".format(volcano_file_path))

    html_file = "index.html"
    html_path = os.path.join(ballgown_output_dir, html_file)
    html_zip_path = html_path + ".zip"
    try:
        f = open(html_path, "w")
    except:
        raise Exception("can't create html file {0}".format(html_path))

    f.write('<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">\n')
    f.write('<html lang="en">\n')
    f.write(' <head>\n')
    f.write('   <meta http-equiv="content-type" content="text/html; charset=utf-8">\n')
    f.write('   <title>title</title>\n')
    f.write(' </head>\n')
    f.write('<body>\n')
    if os.path.exists(volcano_file_path):
        f.write("<img src={0}>\n".format('"' + volcano_plot_file + '"'))
    else:
        f.write('  no volcano plot file generated')
    f.write('  </body>\n')
    f.write('</html>\n')
    f.close()

    with ZipFile(html_zip_path, 'w', allowZip64=True) as izip:
        izip.write(html_path, html_file)

    logger.info("making shock handle for html fiel")
    html_zip_shock_ret = script_util.upload_file_to_shock(logger, html_zip_path)

    logger.info(pformat(html_zip_shock_ret))

    # logger.info( "making shock handle")
    # html_zip_shock_id = script_util.upload_file_to_shock( logger, html_zip_path )['handle']['id']

    # Begin hack.  Giving up on KBaseReport->create_complex_report()
    # Lets roll our own:
    logger.info("hacking KBaseReport.report object")
    kbo = {
        "text_message": "Ballgown FPKM Differential Expression Results",  # string text_message;
        #           "warnings"              : [""],                                         # list<string> warnings;
        "objects_created": [  # list<WorkspaceObject> objects_created;
            {
                "ref": de_obj_ref,
                "description": "Differential expression object"
            },
            {
                "ref": em_obj_ref,
                "description": "Filtered expression Matrix object"
            },
        ],
        # "file_links"            : [                                             # list<LinkedFile> file_links;
        #                            {                                           # LinkedFile;
        #                              #"handle"      : image_zip_shock_ret['handle']['hid'],         # handle_ref handle;   string?
        #                              "handle"      : volcano_file_shock_ret['handle']['hid'],         # handle_ref handle;   string?
        #                              "description" : "volcano plot png file",  # string description;
        #                              "name"        : volcano_plot_file,        # string name;
        #                              #"URL"         : image_zip_shock_ret['handle']['url'] + "/node/" + image_zip_shock_ret['handle']['id']
        #                              "URL"         : volcano_file_shock_ret['handle']['url'] + "/node/" + volcano_file_shock_ret['handle']['id']
        #                            }
        #                          ],
        # "html_links"            : [                                              # list<LinkedFile> html_links;
        #                            {                                            # LinkedFile;
        #                              "handle"      : html_zip_shock_ret['handle']['hid'],          # handle_ref handle;
        #                              "description" : "volcano plot html file",  # string description;
        #                              "name"        : html_file,                 # string name;
        #                              "label"       : None,
        #                              "URL"         : html_zip_shock_ret['handle']['url'] + "/node/" + html_zip_shock_ret['handle']['id']
        #                            }
        #                          ],
        "direct_html": None,  # string direct_html;
        "direct_html_link_index": None  # int direct_html_link_index;
    }
    logger.info(pformat(kbo))
    objs_save_data = ws_client.save_objects(
        {"workspace": ws_id,
         "objects": [
             {
                 "type": "KBaseReport.Report",
                 "data": kbo,
                 "name": report_obj_name,
                 "hidden": 1
             }
         ]
         }
    )
    return objs_save_data[0]

    # logger.info( "initializing report object with callback {0}".format( callback_url ))
    # kbr = KBaseReport( callback_url, token=token)
    # logger.info( pformat( kbr ) )
    # report_input_params = {
    #                        'objects_created'   : [ { 'ref': de_obj_ref, 'description': 'Differential Expression' },
    #                                                { 'ref': em_obj_ref, 'description': 'Filtered Expression Matrix' },
    #                                              ],
    #                        'direct_html_index' : 0,
    #                        'file_links'        : [
    #                                                {
    #                                                  'path'       :  image_zip_path,
    #                                                  'name'       :  image_zip_file,
    #                                                  'description': 'zip file containing volcano plot'
    #                                                }
    #                                              ],
    #                        'html_links'        : [
    #                                                {
    #                                                 'path'       : html_path,
    #                                                 'name'       : html_file,
    #                                                 'description': 'HTML file to display volcano plot'
    #                                                }
    #                                              ],
    #                        'report_object_name': report_obj_name,
    #                        'workspace_name'    : ws_id
    #                      }
    # logger.info( "KBaseReport initialized, trying to upload report, params are" )
    # logger.info( pformat( report_input_params ))
    #
    ##new_input_params = {
    #                     'message' :  "This is a test I hope it works"
    #                   }

    # try:
    #    repout = kbr.create_extended_report( report_input_params )
    #    #repout = kbr.create_extended_report( new_input_params )
    # except:
    #    raise Exception( "Unable to create_extended_report" )
    #
    # return repout['name']
