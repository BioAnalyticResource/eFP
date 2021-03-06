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

        if self.conf['DB_ANNO'] is None or self.conf['DB_LOOKUP_TABLE'] is None or \
                self.conf['LOOKUP'][self.database] == '0':
            # annotations db not defined
            self.retrieve_lookup_gene_data(efp_id)
            return

        if self.conn is None:
            self.connect(self.conf['DB_ANNO'])

        cursor = self.conn.cursor()

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
        return

    def retrieve_lookup_gene_data(self, id):
        """
        Retrieve_lookup_gene_data
        Checks whether a gene exists when no lookup is available e.g. RNA-seq Data
        @param id:
        @return:
        """
        if id is None:
            return

        if self.conn is None:
            self.connect(self.database)

        cursor = self.conn.cursor()

        cursor.execute('SELECT data_probeset_id FROM sample_data WHERE data_probeset_id LIKE %s', (id + '%',))
        row = cursor.fetchall()
        cursor.close()
        self.conn = None

        if len(row) > 0:
            self.gene_id = id
            self.probeset_id = id

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
                     (self.conf['DB_NCBI_GENEID_COL'], self.conf['DB_NCBI_GENE_TABLE'])
        cursor.execute(select_cmd, (ncbi_gi,))
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
            
            # Maize stuff. Might need to figure this out
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
            
            if result:
                self.annotation = result[0]
                cursor.close()

                items = '__'.split(self.annotation)
                aliases = '_'.split(items[0])

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
