import logging
import os
import pandas as pd
import sys
try:
    from utilities.utilities import _TCGA_CANCER_TYPES, _PHEN_TABLE, _PER
except ModuleNotFoundError:
    sys.path.append(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    )
    from utilities.utilities import _TCGA_CANCER_TYPES, _PHEN_TABLE, _PER
from utilities.utilities import _JX_SAMP_TABLE, _CANCER_TYPES_PRIMARY
from utilities.utilities import _JX_ANN_TABLE, _MATCHED_NORMALS


def count_total_jxs(out_path, now, db_conn):
    """

    :param out_path:
    :param now:
    :param db_conn:
    :return:
    """
    logging.info('starting count_total_jxs function')
    sample_count_query = (
        "SELECT samp_count, COUNT (jx_id) jx_count "
        "FROM (SELECT jx_id, COUNT({js}.recount_id) samp_count "
        "FROM {js} INNER JOIN {sp} ON {js}.recount_id == {sp}.recount_id "
        "AND {sp}.tumor_normal==0 GROUP BY jx_id) GROUP BY samp_count;"
        "".format(js=_JX_SAMP_TABLE, sp=_PHEN_TABLE)
    )
    jxs_per_samp_count = pd.read_sql_query(sample_count_query, db_conn)
    file_name = 'jx_counts_per_sample_count_{}.txt'.format(now)
    file_path = os.path.join(out_path, file_name)
    with open(file_path, 'w') as output:
        jxs_per_samp_count.to_csv(output, index=False)

    jx_count_query = (
        'SELECT COUNT (DISTINCT jx_id) tumor_jx_count '
        'FROM {js} INNER JOIN {sp} '
        'ON {js}.recount_id == {sp}.recount_id AND {sp}.tumor_normal == 0;'
        ''.format(js=_JX_SAMP_TABLE, sp=_PHEN_TABLE)
    )
    total_cancer_jxs = pd.read_sql_query(jx_count_query, db_conn)
    total_cancer_jxs = total_cancer_jxs['tumor_jx_count'][0]

    file_name = 'total_tcga_tumor_jxs.txt'
    file_path = os.path.join(out_path, file_name)
    with open(file_path, 'a') as output:
        output.write(
            'total number of junctions occurring in TCGA tumor samples: {}\n'
            ''.format(total_cancer_jxs)
        )
    logging.info('total jx count ending\n')
    return


def count_samples(out_path, now, db_conn, can_set='TCGA'):
    """

    :param out_path:
    :param now:
    :param db_conn:
    :return:
    """
    logging.info('starting count_samples function')

    if can_set == 'TCGA':
        label_col = 'project_type_label'
    elif can_set == 'primary':
        label_col = 'primary_type'

    normal_search_command = (
        'SELECT {lc}, COUNT ({lc}) FROM {sp} WHERE tumor_normal == 1 '
        'GROUP BY {lc};'.format(lc=label_col, sp=_PHEN_TABLE)
    )
    normal_counts = pd.read_sql_query(normal_search_command, db_conn)
    norm_names = {
        'COUNT ({})'.format(label_col): 'normal_samp_count',
        label_col: 'tissue'
    }
    normal_counts.rename(columns=norm_names, inplace=True)

    cancer_search_command = (
        'SELECT {lc}, COUNT ({lc}) FROM {sp} WHERE tumor_normal == 0 '
        'GROUP BY {lc};'.format(lc=label_col, sp=_PHEN_TABLE)
    )

    cancer_counts = pd.read_sql_query(cancer_search_command, db_conn)
    canc_names = {
        'COUNT ({})'.format(label_col): 'tumor_samp_count',
        label_col: 'tissue'
    }
    cancer_counts.rename(columns=canc_names, inplace=True)
    new_df = pd.merge(
        normal_counts, cancer_counts, on='tissue', how='outer'
    ).fillna(0)
    new_df = new_df.astype(
        {'tumor_samp_count': int, 'normal_samp_count': int}
    )
    file_name = 'database_sample_counts_{}.csv'.format(now)
    file_path = os.path.join(out_path, file_name)
    with open(file_path, 'w') as output:
        new_df.to_csv(output, index=False)
    print('ending count samples function\n')
    logging.info('ending count samples function\n')
    return


def collect_all_jxs(batch_num, out_path, now, db_conn):
    """Collects all junctions and their annotations for every cancer type.

    :param batch_num:
    :param out_path:
    :param now:
    :param conn:



    returns None
    """
    print('\nstarting collect all jxs function')
    logging.info('starting collect all jxs function')
    all_cancers = list(
        set(_TCGA_CANCER_TYPES).union(set(_CANCER_TYPES_PRIMARY))
    )
    if batch_num:
        try:
            collect_list = [all_cancers[batch_num - 1]]
        except IndexError:
            return
    else:
        collect_list = all_cancers

    for cancer in collect_list:
        logging.info('collecting all jxs for {}:'.format(cancer))
        print('\ncollecting info for {}:'.format(cancer))

        if cancer in _CANCER_TYPES_PRIMARY:
            cancer_col_name = 'primary_type'
        else:
            cancer_col_name = 'project_type_label'

        sample_count_command = (
            'SELECT COUNT (*) FROM {} WHERE tumor_normal == 0 '
            'AND {} == "{}";'.format(_PHEN_TABLE, cancer_col_name, cancer)
        )
        count = pd.read_sql_query(sample_count_command, db_conn)
        count = count['COUNT (*)'][0]
        logging.info('{} samples'.format(count))
        if count == 0:
            logging.info('no samples, continuing\n')
            continue
        per_col = cancer + _PER
        all_jxs_name = '{}_all_jxs_{}.csv'.format(cancer, now)
        select_command = (
         'SELECT {ja}.jx, {ja}.annotation, COUNT (phen_recount) '
         'FROM (SELECT {js}.jx_id can_jxs, phen_recount '
         'FROM (SELECT recount_id phen_recount '
         'FROM {sp} '
         'WHERE {sp}.{can_col} == "{can}" AND {sp}.tumor_normal == 0) '
         'INNER JOIN {js} ON phen_recount == {js}.recount_id) '
         'INNER JOIN {ja} ON {ja}.jx_id==can_jxs '
         'GROUP BY ({ja}.jx_id);'.format(
             js=_JX_SAMP_TABLE, ja=_JX_ANN_TABLE, sp=_PHEN_TABLE, can=cancer,
             can_col=cancer_col_name
         )
        )
        query_result = pd.read_sql_query(select_command, db_conn)
        col_rename = {'COUNT (phen_recount)': per_col}
        query_result.rename(columns=col_rename, inplace=True)
        query_result[per_col] = query_result[per_col] / count
        query_result = query_result.sort_values(by=[per_col], ascending=False)

        all_jxs_file = os.path.join(out_path, all_jxs_name)
        with open(all_jxs_file, 'w') as output:
            query_result.to_csv(output, index=False)

    print('\nending collect all jxs function\n')
    logging.info('ending collect all jxs function\n')
    return


def create_nonnormal_jx_query(non_gtex, ann_filter, cov_filter, jx_recount_id,
                              cancer):
    """

    :param non_gtex:
    :param ann_filter:
    :param cov_filter:
    :param jx_recount_id:
    :param cancer:
    :return:
    """
    if cancer in _CANCER_TYPES_PRIMARY:
        cancer_col_name = 'primary_type'
    else:
        cancer_col_name = 'project_type_label'

    if jx_recount_id:
        jx_to_select = "jx_id, "
    else:
        jx_to_select = "jx, "

    if ann_filter:
        ann_add = " AND {ja}.annotation < 3".format(ja=_JX_ANN_TABLE)
    else:
        ann_add = ''

    if cov_filter:
        cov_add = " AND cov >= 5"
    else:
        cov_add = ''

    if non_gtex:
        non_normal =  "AND {sp}.tumor_normal == 1) ".format(sp=_PHEN_TABLE)
    else:
        match = _MATCHED_NORMALS[cancer]
        if len(match) > 0:
            normals = '"' + '", "'.join(match) + '"'
            add_match = (
                ' OR ({sp}.{can_col} IN ({nor}))'.format(
                    sp=_PHEN_TABLE, nor=normals, can_col=cancer_col_name
                )
            )
        else:
            add_match = ''
        non_normal = (
            "AND (({sp}.project_type_label == '{can}' "
            "AND {sp}.tumor_normal == 1){ma})) ".format(
                sp=_PHEN_TABLE, can=cancer, ma=add_match
            )
        )

    count_query = (
        "SELECT DISTINCT {jta}jx_recounts recount_id, tcga_id, norm_jxs "
        "  FROM (SELECT jx, jx_id, jx_recounts, {sp}.tcga_id "
        "    FROM (SELECT jx, {js}.jx_id, {js}.recount_id jx_recounts, "
        "          {js}.coverage cov "
        "          FROM {js} INNER JOIN {ja} ON {js}.jx_id == {ja}.jx_id{aa}) "
        "    INNER JOIN {sp} ON jx_recounts == {sp}.recount_id "
        "      AND ({sp}.{can_col} == '{can}'{ca})) "
        "  LEFT JOIN (SELECT DISTINCT {js}.jx_id norm_jxs FROM {js} "
        "    INNER JOIN {sp} ON {js}.recount_id == {sp}.recount_id {nn}"
        "ON jx_id == norm_jxs WHERE norm_jxs IS NULL;".format(
            jta=jx_to_select, js=_JX_SAMP_TABLE, aa=ann_add, ca=cov_add,
            sp=_PHEN_TABLE, can_col=cancer_col_name, nn=non_normal, can=cancer,
            ja=_JX_ANN_TABLE
        )
    )
    return count_query


def collect_sample_jxs(cancer, out_path, now, db_conn, normal=True):
    """

    :param cancer:
    :param out_path:
    :param now:
    :param db_conn:
    :param normal:
    :return:
    """
    if normal:
        t_flag = 'normal'
    else:
        t_flag = 'tumor'

    all_jxs_name = (
        '{}_all_jxs_{}_NOcov_NOann_filter_{}.txt'.format(cancer, t_flag, now)
    )

    if cancer in _CANCER_TYPES_PRIMARY:
        cancer_col_name = 'primary_type'
    else:
        cancer_col_name = 'project_type_label'

    sql_jx_query = (
        "SELECT DISTINCT jx, jx_recounts recount_id, tcga_id "
        "  FROM (SELECT jx, jx_id, jx_recounts, {sp}.tcga_id "
        "    FROM (SELECT jx, {js}.jx_id, {js}.recount_id jx_recounts, "
        "          {js}.coverage cov "
        "          FROM {js} INNER JOIN {ja} ON {js}.jx_id == {ja}.jx_id) "
        "    INNER JOIN {sp} ON jx_recounts == {sp}.recount_id "
        "      AND {sp}.{can_col} == '{can}'"
        "      AND {sp}.tumor_normal == {tf}"
        "      );".format(
            js=_JX_SAMP_TABLE, sp=_PHEN_TABLE, can_col=cancer_col_name,
            can=cancer, ja=_JX_ANN_TABLE, tf=int(normal)
        )
    )

    logging.info('timestamp is: {}'.format(now))
    logging.info(sql_jx_query)
    neojx_count = pd.read_sql_query(sql_jx_query, db_conn)

    all_jxs_file = os.path.join(out_path, all_jxs_name)
    with open(all_jxs_file, 'w') as output:
        neojx_count.to_csv(output, index=False)

    return


def count_and_collect_neojxs(batch_num, out_path, now, db_conn, non_gtex=False,
                             cov_filter=True, ann_filter=False, print_jxs=True,
                             print_counts=True, addl_flag='',
                             jx_recount_id=True):
    """

    Input:
    jx_recount_id: whether to output only the recount ID ("true") or also the
        TCGA ID ("False") (Boolean)
    two_filters: whether to filter with minimum 5 scaled jx coverage (True) or
        no coverage filter (False) (Boolean)
    non_gtex: whether to filter out all GTEx junctions (True) or only junctions
        from non-paired normal GTEx tissue (False) (Boolean)

    :param batch_num:
    :param out_path:
    :param now:
    :param db_conn:
    :param non_gtex:
    :return:
    """
    logging.info('starting neojxs both function')

    if non_gtex:
        normals = 'GTEx'
        out_flag = 'non_GTEx'
        dir_flag = 'non-core-normal'
    else:
        out_flag = 'non_paired_normal'
        dir_flag = 'non-tissue-matched'

    if cov_filter:
        cov_flag = ''
    else:
        cov_flag = 'NO'

    if ann_filter:
        ann_flag = ''
    else:
        ann_flag = 'NO'

    if jx_recount_id:
        jx_flag = ''
    else:
        jx_flag = 'jx_'

    logging.info('counting total number of unique neojxs in TCGA:')
    if ann_filter:
        ann_addition = ' AND {ja}.annotation < 3'
    else:
        ann_addition = ''

    count_unique_query = (
        "SELECT COUNT (DISTINCT can_jx) neojx_count FROM "
        "    (SELECT can_jx FROM (SELECT DISTINCT {js}.jx_id can_jx "
        "        FROM {js} INNER JOIN {ja} "
        "        ON {js}.jx_id == {ja}.jx_id{ann}) "
        "    LEFT JOIN "
        "        (SELECT DISTINCT {js}.jx_id norm_jxs "
        "        FROM {js} INNER JOIN {sp} "
        "        ON {js}.recount_id == {sp}.recount_id "
        "        AND {sp}.tumor_normal == 1) "
        "    ON can_jx == norm_jxs WHERE norm_jxs IS NULL);".format(
            js=_JX_SAMP_TABLE, sp=_PHEN_TABLE, ja=_JX_ANN_TABLE,
            ann=ann_addition
        )
    )

    unique_neojx_count = pd.read_sql_query(count_unique_query, db_conn)
    logging.info(
        'total number of unique neojxs in TCGA: {}'
        ''.format(unique_neojx_count['neojx_count'])
    )
    logging.info('collecing counts per sample:')

    count_name = 'tcga_total_neojx_counts_'
    if batch_num:
        count_name += 'batch'
        full_list = list(
            set(_TCGA_CANCER_TYPES).union(set(_CANCER_TYPES_PRIMARY))
        )
        try:
            collect_list = [full_list[batch_num - 1]]
        except IndexError:
            logging.info('collect list error:')
            logging.info(full_list)
            logging.info(len(full_list))
            return
    else:
        collect_list = _MATCHED_NORMALS

    count_name += '{}_{}.txt'.format(now, out_flag)
    count_file = os.path.join(out_path, count_name)

    all_neojxs = set()
    for cancer in collect_list:
        logging.info('collecting neojxs for {}:'.format(cancer))

        if cancer in _CANCER_TYPES_PRIMARY:
            cancer_col_name = 'primary_type'
        else:
            cancer_col_name = 'project_type_label'

        sample_count_command = (
            'SELECT COUNT (*) FROM {} WHERE tumor_normal == 0 '
            'AND {} == "{}";'.format(_PHEN_TABLE, cancer_col_name, cancer)
        )
        count = pd.read_sql_query(sample_count_command, db_conn)
        count = count['COUNT (*)'][0]
        logging.info('{} samples'.format(count))
        if count == 0:
            logging.info('no samples, continuing')
            continue
        logging.info('collecting info for {}:'.format(cancer))

        if not non_gtex:
            try:
                match = _MATCHED_NORMALS[cancer]
            except KeyError:
                logging.info('no matched normals; non-TCGA cancer. continuing')
                continue
            normals = ', '.join(match)

        jx_per_samp_name = (
            '{}{}_neojx_counts_per_sample_{}_{}_{}cov_{}ann_filter.txt'
            ''.format(jx_flag, cancer, now, out_flag, cov_flag, ann_flag)
        )
        all_neojxs_name = (
            '{}{}_all_neojxs_{}_{}_{}cov_{}ann_filter.txt'
            ''.format(jx_flag, cancer, now, out_flag, cov_flag, ann_flag)
        )

        sql_count_query = create_nonnormal_jx_query(
            non_gtex, ann_filter, cov_filter, jx_recount_id, cancer
        )
        logging.info('timestamp is: {}'.format(now))
        logging.info(sql_count_query)
        print(sql_count_query)
        neojx_count = pd.read_sql_query(sql_count_query, db_conn)
        print(neojx_count, '\n')
        neojx_count.drop(['norm_jxs'], inplace=True, axis=1)

        if print_jxs:
            all_jxs_dir = os.path.join(
                out_path, '{}{}_all_jxs_per_sample'.format(addl_flag, dir_flag)
            )
            os.makedirs(all_jxs_dir, exist_ok=True)
            all_neojxs_file = os.path.join(all_jxs_dir, all_neojxs_name)
            with open(all_neojxs_file, 'w') as output:
                neojx_count.to_csv(output, index=False)

        if jx_recount_id:
            unique_jxs = neojx_count.jx_id.unique()
            tot_neojxs = len(neojx_count.jx_id)

        else:
            unique_jxs = neojx_count.jx.unique()
            tot_neojxs = len(neojx_count.jx)
        unique_neojxs = len(unique_jxs)

        if not batch_num:
            all_neojxs = all_neojxs.union(set(unique_jxs.tolist()))

        if jx_recount_id:
            tempgroup = neojx_count.groupby(['recount_id'])['jx_id'].count()
        else:
            tempgroup = neojx_count.groupby(['recount_id'])['jx'].count()
        grouped_neojxs = pd.DataFrame(
            {'recount_id': tempgroup.index,
             'neojx_count': tempgroup.values}
        )
        mean_count = grouped_neojxs.neojx_count.mean()

        with open(count_file, 'a') as output:
            output.write(
                '{} not in {}: {} total jxs, {} unique jxs, {} samples, '
                'avg {} jxs/sample\n'
                ''.format(
                    cancer, normals, tot_neojxs, unique_neojxs, count,
                    mean_count
                )
            )

        if print_counts:
            count_dir = os.path.join(
                out_path, '{}{}_counts_per_sample'.format(addl_flag, dir_flag)
            )
            os.makedirs(count_dir, exist_ok=True)
            jx_per_samp_file = os.path.join(count_dir, jx_per_samp_name)
            with open(jx_per_samp_file, 'w') as output:
                grouped_neojxs.to_csv(output, index=False)

    if not batch_num:
        jx_count_query = (
            'SELECT COUNT (DISTINCT jx_id) tumor_jx_count '
            'FROM {js} INNER JOIN {sp} '
            '   ON {js}.recount_id == {sp}.recount_id '
            '   AND {sp}.tumor_normal == 0;'
            ''.format(js=_JX_SAMP_TABLE, sp=_PHEN_TABLE)
        )
        total_cancer_jxs = pd.read_sql_query(jx_count_query, db_conn)
        total_cancer_jxs = total_cancer_jxs['tumor_jx_count'][0]
        with open(count_file, 'a') as output:
            output.write(
                '\ntotal number of junctions occurring '
                'in all TCGA tumor samples: {}'
                ''.format(total_cancer_jxs)
            )
            output.write(
                '\ntotal unique neojxs in TCGA: {}\n'.format(len(all_neojxs))
            )
    print('\nending count neojxs both function\n')
    logging.info('ending count neojxs both function\n')
    return


def non_normal_jxs_prevs(batch_num, out_path, now, db_conn, index_db,
                         tcga_types=False, all_types=False, ann_filter=False):
    """

    :param batch_num:
    :param out_path:
    :param now:
    :param db_conn:
    :param index_db:
    :param tcga_types:
    :param all_types:
    :param ann_filter:
    :return:
    """
    print('\nstarting non-normal jx collection')

    norm_temp_command = (
        'CREATE TEMP TABLE normal AS SELECT DISTINCT jx_id FROM {} '
        'INNER JOIN {} ON {}.recount_id == {}.recount_id '
        'WHERE tumor_normal == 1;'
        ''.format(_JX_SAMP_TABLE, _PHEN_TABLE, _JX_SAMP_TABLE, _PHEN_TABLE)
    )
    index_db.execute(norm_temp_command)
    db_conn.commit()

    if tcga_types:
        cancer_list = _TCGA_CANCER_TYPES
    elif all_types:
        cancer_list = list(
            set(_TCGA_CANCER_TYPES).union(set(_CANCER_TYPES_PRIMARY))
        )
    else:
        cancer_list = _CANCER_TYPES_PRIMARY

    if batch_num:
        collect_list = [cancer_list[batch_num - 1]]
    else:
        collect_list = cancer_list

    for cancer in collect_list:
        if cancer in _CANCER_TYPES_PRIMARY:
            cancer_col_name = 'primary_type'
        else:
            cancer_col_name = 'project_type_label'
        per_col = cancer + _PER
        sample_count_command = (
            'SELECT COUNT (*) FROM {} WHERE tumor_normal == 0 '
            'AND {} == "{}";'.format(_PHEN_TABLE, cancer_col_name, cancer)
        )
        count = pd.read_sql_query(sample_count_command, db_conn)
        count = count['COUNT (*)'][0]
        if count == 0:
            continue
        print('\ncollecting info for {}:'.format(cancer))

        if ann_filter:
            ann_addition = ' AND {ja}.annotation < 3'
            ann_flag = 'unannotated_'
        else:
            ann_addition = ''
            ann_flag = ''

        select_command = (
            'SELECT {ja}.jx, {ja}.annotation, COUNT (phen_recount) '
            'FROM (SELECT phen_recount, nonnorm_jxs '
            'FROM (SELECT {js}.jx_id nonnorm_jxs, phen_recount '
            'FROM (SELECT recount_id phen_recount '
            'FROM {sp} '
            'WHERE {sp}.{can_col} == "{can}" AND {sp}.tumor_normal == 0) '
            'INNER JOIN {js} ON phen_recount == {js}.recount_id) '
            'LEFT OUTER JOIN temp.normal nor ON nonnorm_jxs == nor.jx_id '
            'WHERE nor.jx_id IS NULL) '
            'INNER JOIN {ja} '
            'ON {ja}.jx_id==nonnorm_jxs{ann_add} '
            'GROUP BY ({ja}.jx_id);'.format(
                js=_JX_SAMP_TABLE, ja=_JX_ANN_TABLE, sp=_PHEN_TABLE,
                can=cancer, can_col=cancer_col_name, ann_add=ann_addition
            )
        )

        query_result = pd.read_sql_query(select_command, db_conn)
        col_rename = {'COUNT (phen_recount)': per_col}
        query_result.rename(columns=col_rename, inplace=True)
        query_result[per_col] = query_result[per_col] / count
        query_result = query_result.sort_values(by=[per_col], ascending=False)

        jx_filename = (
            '{}_not-in-GTEx_{}jxs_{}.csv'.format(cancer, ann_flag, now)
        )
        full_path = os.path.join(out_path, jx_filename)
        with open(full_path, 'w') as output:
            query_result.to_csv(output, index=False)

    index_db.execute('DROP TABLE IF EXISTS temp.normal;')
    db_conn.commit()
    return


def collect_data_for_analyses(batch_num, out_path, now, conn, index_db):
    """Completes all queries and file generation for neojx paper analyses

    :param batch_num:
    :param out_path:
    :param now:
    :param conn:
    :param index_db:
    :return:
    """
    # Count samples for Table S1
    count_samples(out_path, now, conn, can_set='primary')
    count_samples(out_path, now, conn)

    # Collect all junctions for all TCGA samples for set membership annotation
    all_jxs_dir = os.path.join(out_path, 'all_jxs')
    os.makedirs(all_jxs_dir, exist_ok=True)
    collect_all_jxs(batch_num, all_jxs_dir, now, conn)

    # For use in Figs S3 and S4 and set membership annotation
    ncn_dir = os.path.join(out_path, 'non-core-normal_jxs_per_sample')
    os.makedirs(ncn_dir, exist_ok=True)
    count_and_collect_neojxs(
        batch_num, ncn_dir, now, conn, non_gtex=True, cov_filter=False,
        ann_filter=False, jx_recount_id=False
    )
    # For use in set membership annotation
    ntm_dir = os.path.join(out_path, 'non-tissue-matched_jxs_per_sample')
    os.makedirs(ntm_dir, exist_ok=True)
    count_and_collect_neojxs(
        batch_num, ntm_dir, now, conn, non_gtex=False, cov_filter=False,
        ann_filter=False, jx_recount_id=False, print_counts=False
    )

    # For use in Figures 2A, 2B, 2C, and S5
    prev_dir = os.path.join(out_path, 'non-core-normal_jx_prevalences')
    os.makedirs(prev_dir, exist_ok=True)
    non_normal_jxs_prevs(
        batch_num, prev_dir, now, conn, index_db, all_types=True
    )

    # For use in Figure S7
    norm_dir = os.path.join(out_path, 'all_jxs_per_sample_paired_normals')
    os.makedirs(norm_dir, exist_ok=True)
    collect_sample_jxs(
        'Skin_Cutaneous_Melanoma', norm_dir, now, conn, normal=False
    )
    collect_sample_jxs(
        'Skin', norm_dir, now, conn, normal=True
    )
    collect_sample_jxs(
        'Skin_Cutaneous_Melanoma', norm_dir, now, conn, normal=True
    )

    # For use in Fig S1
    count_and_collect_neojxs(
        batch_num, out_path, now, conn, non_gtex=True, cov_filter=False,
        ann_filter=False, jx_recount_id=False, addl_flag='FigS1_',
        print_jxs=False
    )
    count_and_collect_neojxs(
        batch_num, out_path, now, conn, non_gtex=False, cov_filter=False,
        ann_filter=False, jx_recount_id=False, addl_flag='FigS1_',
        print_jxs=False
    )

    return


def main(args, now, conn, index_db):
    out_path = args.output_path
    batch_num = args.batch_number
    log_mode = args.log_level

    if batch_num:
        log_file = os.path.join(
            out_path, 'query_log_{}.txt'.format(batch_num)
        )
    else:
        log_file = os.path.join(
            out_path, 'query_log_{}.txt'.format(now)
        )
    logging.basicConfig(filename=log_file, level=log_mode)

    collect_data_for_analyses(batch_num, out_path, now, conn, index_db)
