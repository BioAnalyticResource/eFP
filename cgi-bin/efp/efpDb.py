#!/usr/bin/python3
import re
import MySQLdb
import sys


class Gene:
    def __init__(self, gene_id, database, conf):
        self.conn = None
        self.annotation = None
        self.alias = None
        self.gene_id = None
        self.probeset_id = None
        self.ncbi_id = None
        self.database = database
        self.lookup = None
        self.webservice_gene = None
        self.conf = conf

        # eFP Tomato
        if self.conf['species'] == "TOMATO":
            # Only keep isoforms in the following
            if (self.database not in ["tomato_ils", "tomato_ils2", "tomato_s_pennellii", "tomato_seed"]) and (gene_id is not None):
                gene_id = re.sub(r"\.\d$", "", gene_id)

        # eFP Maize
        if self.conf['species'] == "MAIZE":
            # Convert ID to v3 if possible, otherwise return id as is (v4)
            if re.search(r'^Zm\d+d\d+$', gene_id, re.I):
                gene_id = self.maize_v4_to_v3(gene_id)

        self.retrieve_gene_data(gene_id)

        if self.gene_id is None:
            self.ncbi_to_gene_id(gene_id)
            self.retrieve_gene_data(self.gene_id)

    def get_gene_id(self):
        return self.gene_id

    def get_probeset_id(self):
        return self.probeset_id

    def get_ncbi_id(self):
        return self.ncbi_id

    def retrieve_gene_data(self, efp_id):
        """
        retrieve_gene_data: Retrieves the probeset ID that corresponds to the given gene ID
        @param efp_id:
        @return:
        """

        if efp_id is None:
            return

        # This part is only for eFP Maize.
        if self.conf['species'] == "MAIZE":
            if self.database == "maize_gdowns":
                if self.conn is None:
                    self.connect(self.conf['DB_GDOWNS_LOOKUP'])

                cursor = self.conn.cursor()
                select_cmd = "SELECT gene, probeset FROM %s WHERE %s=%%s AND date = (SELECT MAX(date) FROM %s)" % (
                    self.conf['DB_GDOWNS_LOOKUP_TABLE'], self.conf['DB_ANNO_GENEID_COL'], self.conf['DB_GDOWNS_LOOKUP_TABLE'])

                cursor.execute(select_cmd, (efp_id,))
                row = cursor.fetchall()
                cursor.close()
                self.conn = None

                if len(row) > 0:
                    self.gene_id = row[0][0]
                    self.probeset_id = row[0][1]
                    return

            if self.database in ["maize_buell_lab", "maize_early_seed", "maize_embryonic_leaf_development"]:
                self.conn = None
                self.connect(self.conf['DB_ANNO'])
                cursor = self.conn.cursor()

                select_cmd = "SELECT v4, v3 FROM %s WHERE v3=%%s" % (self.conf['DB_BUELL_LOOKUP_TABLE'])

                cursor.execute(select_cmd, (efp_id,))
                row = cursor.fetchall()
                cursor.close()
                self.conn = None
                if len(row) > 0:
                    # Both v4 and v3 exists
                    self.gene_id = row[0][0]
                    self.probeset_id = row[0][0]
                    return
                else:
                    # You are here if v3 version of the v4 gene doesn't exist.
                    # Make sure v4 gene actually exists, and not just some made up thing that matches the format!
                    if self.database == "maize_buell_lab":
                        self.connect("maize_buell_lab")
                    elif self.database == "maize_early_seed":
                        self.connect("maize_early_seed")
                    elif self.database == "maize_embryonic_leaf_development":
                        self.connect("maize_embryonic_leaf_development")
                    cursor = self.conn.cursor()

                    select_cmd = "SELECT data_probeset_id FROM sample_data WHERE data_probeset_id=%s"
                    cursor.execute(select_cmd, (efp_id,))
                    row = cursor.fetchall()
                    cursor.close()
                    self.conn = None

                    if len(row) > 0:
                        self.gene_id = row[0][0]
                        self.probeset_id = None
                        return
                    else:
                        return

        if self.conf['DB_ANNO'] is None or self.conf['DB_LOOKUP_TABLE'] is None or \
                self.conf['LOOKUP'][self.database] == '0':
            # annotations db not defined
            self.retrieve_lookup_gene_data(efp_id)
            return

        if self.conn is None:
            self.connect(self.conf['DB_ANNO'])

        cursor = self.conn.cursor()

        # Special case for Medicago Seed
        if (self.conf['species'] == "MEDICAGO") and (self.database == "medicago_seed"):
            self.conf['DB_LOOKUP_TABLE'] = self.conf["DB_LOOKUP_TABLE_SEED"]

        # First check if the species has an orthologous lookup
        if 'LOOKUP_ORTHO' in self.conf and self.conf['LOOKUP_ORTHO']:
            select_cmd = 'SELECT t1.%s, t1.%s FROM %s t1 WHERE (t1.%s=%%s or t1.%s=%%s) AND t1.date=(SELECT MAX(t2.date) FROM %s t2)' % \
                         (self.conf['DB_LOOKUP_ARABIDOPSIS_COL'], self.conf['DB_LOOKUP_GENEID_COL'],
                          self.conf['DB_LOOKUP_TABLE'],
                          self.conf['DB_LOOKUP_ARABIDOPSIS_COL'], self.conf['DB_LOOKUP_GENEID_COL'],
                          self.conf['DB_LOOKUP_TABLE'])
        elif self.conf['species'] in ["TOMATO"]:
            select_cmd = 'SELECT t1.`%s`, t1.`%s` FROM `%s` t1 WHERE (t1.`%s`=%%s or t1.`%s`=%%s) AND t1.date=(SELECT MAX(t2.date) FROM `%s` t2)' % \
                 (self.conf['DB_LOOKUP_GENEID_COL'], self.conf['DB_LOOKUP_PBSET_COL'], self.conf['DB_LOOKUP_TABLE'],
                  self.conf['DB_LOOKUP_GENEID_COL'], self.conf['DB_LOOKUP_PBSET_COL'], self.conf['DB_LOOKUP_TABLE'])
        elif self.conf['species'] in ["MOUSE"]:
            select_cmd = 'SELECT t1.`%s`, t1.probe_id FROM `%s` t1 WHERE (t1.probe_id=%%s or t1.`%s`=%%s)' % \
                         (self.conf['DB_LOOKUP_GENEID_COL'], self.conf['DB_LOOKUP_TABLE'],
                          self.conf['DB_LOOKUP_GENEID_COL'])
        else:
            select_cmd = 'SELECT t1.`%s`, t1.probeset FROM `%s` t1 WHERE (t1.probeset=%%s or t1.`%s`=%%s) AND t1.date=(SELECT MAX(t2.date) FROM `%s` t2)' % \
                         (self.conf['DB_LOOKUP_GENEID_COL'], self.conf['DB_LOOKUP_TABLE'],
                          self.conf['DB_LOOKUP_GENEID_COL'], self.conf['DB_LOOKUP_TABLE'])

        cursor.execute(select_cmd, (efp_id, efp_id))
        row = cursor.fetchall()
        cursor.close()
        self.conn = None

        if len(row) > 0:
            self.gene_id = row[0][0]
            self.probeset_id = row[0][1]
        else:
            # Check the sample data table as a last resort
            self.retrieve_lookup_gene_data(efp_id)
        return

    def retrieve_lookup_gene_data(self, efp_id):
        """
        Retrieve_lookup_gene_data
        Checks whether a gene exists when no lookup is available e.g. RNA-seq Data
        @param efp_id:
        @return:
        """
        if efp_id is None:
            return

        # This part is only for eFP Maize.
        if self.conf['species'] == "MAIZE":
            if self.database in ["maize_gdowns", "maize_buell_lab", "maize_early_seed", "maize_embryonic_leaf_development"]:
                return

        if self.conn is None:
            self.connect(self.database)

        cursor = self.conn.cursor()

        cursor.execute('SELECT data_probeset_id FROM sample_data WHERE data_probeset_id LIKE %s', (efp_id + '%',))
        row = cursor.fetchall()
        cursor.close()
        self.conn = None

        if len(row) > 0:
            self.gene_id = efp_id
            self.probeset_id = row[0][0]

        return

    def ncbi_to_gene_id(self, ncbi_gi):
        """
        Returns the AGI corresponding to the NCBI gi accession
        NCBI gi accession comes from NCBI Link out. Need to check whether NCBI gi accession is a NCBI GeneID or NCBI RefSeq.
        @param ncbi_gi:
        @return:
        """
        if ncbi_gi is None:
            return None

        if self.conf['DB_ANNO'] is None or self.conf['DB_NCBI_GENE_TABLE'] is None:  # ncbi lookup db not defined
            return None

        if self.conn is None:
            self.connect(self.conf['DB_ANNO'])

        cursor = self.conn.cursor()

        select_cmd = 'SELECT t1.`%s` FROM `%s` t1 WHERE t1.geneid=%%s or t1.protid=%%s' % \
                     (self.conf['DB_NCBI_GENEID_COL'], self.conf['DB_NCBI_ID_TABLE'])
        cursor.execute(select_cmd, (ncbi_gi, ncbi_gi))
        row = cursor.fetchall()
        cursor.close()

        if len(row) != 0:
            self.gene_id = row[0][0]
            self.ncbi_id = ncbi_gi
        return

    def get_lookup(self):
        if self.database == "maize_rice_comparison":
            if self.conn is None:
                self.connect(self.conf.DB_ANNO)
            cursor = self.conn.cursor()

            maize_convert3 = re.match(r"GRMZM2G[0-9]{6}_T[0-9]{1,2}", self.gene_id)
            if maize_convert3 is not None:
                self.gene_id = re.sub("_T[0-9]{1,2}", "", self.gene_id)

            maize_convert4 = re.match(r"^AC[0-9]{6}\.[0-9]_FGT[0-9]{3}$", self.gene_id)
            if maize_convert4 is not None:
                self.gene_id = self.gene_id.replace("FGT", "FG")

            select_cmd = 'SELECT rice_id FROM `%s` WHERE `%s`=%%s AND date = (SELECT MAX(date) FROM `%s`)' % \
                         (self.conf['DB_ORTHO_LOOKUP_TABLE'], self.conf['DB_ORTHO_GENEID_COL'],
                          self.conf['DB_ORTHO_LOOKUP_TABLE'])
            cursor.execute(select_cmd, (self.gene_id,))
            result = cursor.fetchone()

            if result:
                self.lookup = result[0]
                cursor.close()
            return self.lookup

        if self.database == "rice_maize_comparison":
            if self.conn is None:
                self.connect(self.conf.DB_ANNO)
            cursor = self.conn.cursor()

            select_cmd = "SELECT maize_id FROM `%s` WHERE `%s`=%%s AND date = (SELECT MAX(date) FROM `%s`)" % \
                         (self.conf['DB_ORTHO_LOOKUP_TABLE'], self.conf['DB_ORTHO_GENEID_COL'],
                          self.conf['DB_ORTHO_LOOKUP_TABLE'])
            cursor.execute(select_cmd, (self.gene_id,))
            result = cursor.fetchone()

            if result:
                self.lookup = result[0]
                cursor.close()
            return self.lookup

    def get_annotation(self):
        if self.conf['DB_ANNO'] is None or self.conf['DB_ANNO_TABLE'] is None:
            # annotations db not defined
            return None

        if self.annotation is None:
            if self.conn is None:
                self.connect(self.conf['DB_ANNO'])

            cursor = self.conn.cursor()

            if self.conf['species'] == 'MAIZE':
                maize_convert1 = re.match(r"^GRMZM2G[0-9]{6}$", self.gene_id)
                maize_convert2 = re.match(r"^AC[0-9]{6}\.[0-9]_FG[0-9]{3}$", self.gene_id)
                maize_sp1 = "_T01"

                if maize_convert1 is not None:
                    self.gene_id = self.gene_id + maize_sp1

                if maize_convert2 is not None:
                    self.gene_id = self.gene_id.replace("FG", "FGT")

            select_cmd = "SELECT annotation FROM `%s` WHERE `%s`=%%s AND date = (SELECT MAX(date) FROM `%s`)" % \
                         (self.conf['DB_ANNO_TABLE'], self.conf['DB_ANNO_GENEID_COL'], self.conf['DB_ANNO_TABLE'])

            # Exception for Human
            if self.conf['species'] == 'HUMAN':
                cursor.execute(select_cmd, (self.probeset_id,))
            else:
                cursor.execute(select_cmd, (self.gene_id,))

            result = cursor.fetchone()

            # Initialize
            items = ['']
            aliases = ['']

            if result:
                self.annotation = result[0]
                cursor.close()

                try:
                    items = '__'.split(self.annotation)
                    aliases = '_'.split(items[0])
                except ValueError:
                    if self.annotation == '' or self.annotation is None:
                        items = ['']
                        aliases = ['']

                if len(items) == 1:
                    aliases[0] = ''
                self.alias = aliases[0]

        return self.annotation

    def get_alias(self):
        if self.alias is None:
            self.get_annotation()
        return self.alias

    def connect(self, db_name):
        try:
            self.conn = MySQLdb.connect(host=self.conf['DB_HOST'], user=self.conf['DB_USER'],
                                        passwd=self.conf['DB_PASSWD'], db=db_name)
        except MySQLdb.Error as e:
            print("Error %d: %s" % (e.args[0], e.args[1]), file=sys.stderr)

    def get_homoeologues(self):
        """
        This is for efp Wheat only.
        """

        # This function uses a species connection to wheat database
        try:
            wheat_conn = MySQLdb.connect(host=self.conf['DB_HOST'], user=self.conf['DB_USER'],
                                        passwd=self.conf['DB_PASSWD'], db="wheat")
        except MySQLdb.Error as e:
            print("Error %d: %s" % (e.args[0], e.args[1]), file=sys.stderr)
            exit()

        cursor = wheat_conn.cursor()

        homeologs = list()
        sql = "SELECT gene_id, homeologues from homeologues_lookup where gene_id = %s"

        cursor.execute(sql, (self.gene_id,))
        row = cursor.fetchall()
        cursor.close()
        wheat_conn.close()

        if len(row) > 0:
            if row[0][1] is not None:
                homeologs.extend(row[0][1].split(","))

        return homeologs

    def maize_v4_to_v3(self, gene_id):
        """
        This converts maize v4 gene to v3.
        This is for eFP Maize only.

        @param gene_id:
        @return:
        """
        # This function uses a species connection to wheat database
        try:
            maize_conn = MySQLdb.connect(host=self.conf['DB_HOST'], user=self.conf['DB_USER'],
                                        passwd=self.conf['DB_PASSWD'], db=self.conf['DB_ANNO'])
        except MySQLdb.Error as e:
            print("Error %d: %s" % (e.args[0], e.args[1]), file=sys.stderr)
            exit()

        cursor = maize_conn.cursor()

        select_cmd = "SELECT v3 FROM %s WHERE v4=%%s" % (self.conf['DB_BUELL_LOOKUP_TABLE'])
        cursor.execute(select_cmd, (gene_id,))
        row = cursor.fetchall()
        cursor.close()
        maize_conn.close()

        # The length of row will most likely be > 0
        # But if the v3 column is null, return the original id.
        if len(row) > 0:
            if row[0][0] is not None:
                return row[0][0]
            else:
                return gene_id
        else:
            return gene_id

