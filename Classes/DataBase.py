import sqlite3
from Classes.RDkit import mol_mass_from_smiles, iupac_from_smiles, return_morganfp, similiaryty_list_return


class DataBase:
    def __init__(self, path='./external/resources/data.db'):
        self.conn = sqlite3.connect(path)
        self.cur = self.conn.cursor()
        self.create_table()

    def create_table(self):
        self.cur.execute(
            '''CREATE TABLE IF NOT EXISTS compounds(          
                id INTEGER PRIMARY KEY,
                name TEXT NOT NULL,
                smiles TEXT NOT NULL,
                solutor_id integer,
                molecular_weight REAL,
                iupac_name text,
                morganfp blob,
                FOREIGN KEY (solutor_id) REFERENCES solutors(id)
                )''')

        self.cur.execute(
            '''CREATE TABLE IF NOT EXISTS solutors(          
            id INTEGER PRIMARY KEY,
            name TEXT NOT NULL UNIQUE, 
            UNIQUE ("name") ON CONFLICT IGNORE
            )''')

        self.cur.execute(
            '''CREATE TABLE IF NOT EXISTS cell_lines(          
            id INTEGER PRIMARY KEY,
            name TEXT NOT NULL UNIQUE, 
            description Text,
            UNIQUE ("name") ON CONFLICT IGNORE
            )''')

        self.cur.execute(
            '''CREATE TABLE IF NOT EXISTS MTT(          
            id INTEGER PRIMARY KEY,
            date text,
            compound_id integer REFERENCES compounds (id),
            cell_line_id integer REFERENCES cell_lines (id),
            EC50 real ,
            SE real
            )''')



        self.conn.commit()

    def add_compound(self, name, smiles, morganfp=None, solutor_id=1, molecular_weight=0.0):
        try:
            molecular_weight = mol_mass_from_smiles(smiles)
        except:
            pass

        try:
            iupac_name = iupac_from_smiles(smiles)
        except:
            iupac_name = ''

        self.cur.execute("INSERT INTO compounds(name, smiles, solutor_id, molecular_weight, iupac_name, morganfp) "
                         "VALUES(?,?,?,?,?,?)",
                         (name, smiles,solutor_id, molecular_weight, iupac_name, morganfp))
        self.conn.commit()
        return self.cur.lastrowid

    def add_solutor(self, name):
        self.cur.execute(f"INSERT INTO solutors(name) VALUES(?)",
                         (name,))
        self.conn.commit()
        return self.cur.lastrowid

    def add_cellLine(self, name, description):
        self.cur.execute(f"INSERT INTO cell_lines(name, description) VALUES(?,?)",
                         (name, description, ))
        self.conn.commit()
        return self.cur.lastrowid

    def add_MTT(self, date, compound_id, cell_line_id, EC50, SE):
        self.cur.execute(f"INSERT INTO MTT(date, compound_id, cell_line_id, EC50, SE) VALUES(?,?,?,?,?)",
                         (date, compound_id, cell_line_id, EC50, SE))
        self.conn.commit()
        return self.cur.lastrowid

    def show_compounds(self):
        return self.cur.execute('''SELECT * FROM compounds 
        LEFT JOIN solutors ON compounds.solutor_id = solutors.id 
        ORDER BY id ASC''').fetchall()

    def show_smile_by_id(self, id):
        return self.cur.execute('''SELECT smiles FROM compounds 
        where id = ?''', (id,)).fetchone()

    def show_name_by_id(self, id):
        return self.cur.execute('''SELECT name FROM compounds 
        where id = ?''', (id,)).fetchone()

    def show_solutors(self):
        return self.cur.execute('''SELECT * FROM solutors 
        ORDER BY id ASC''').fetchall()

    def show_mtt(self, comp_id=None):
        if comp_id == None:
            return self.cur.execute('''SELECT * FROM MTT 
                LEFT JOIN compounds ON MTT.compound_id = compounds.id
                LEFT JOIN cell_lines ON MTT.cell_line_id = cell_lines.id
                ''').fetchall()

        else:
            return self.cur.execute('''SELECT * FROM MTT 
                LEFT JOIN compounds ON MTT.compound_id = compounds.id
                LEFT JOIN cell_lines ON MTT.cell_line_id = cell_lines.id
                where compound_id = ?''', (comp_id,)).fetchall()


    def show_lines(self):
        return self.cur.execute('''SELECT * FROM cell_lines 
        ORDER BY id ASC''').fetchall()

    # def show_compound_by_id(self, id):
    #     return self.cur.execute('''SELECT * FROM compounds ORDER BY id ASC''').fetchall()

    def del_compound_by_id(self, id):
        self.cur.execute(f'DELETE from compounds where id = {id}')
        self.conn.commit()

    def showMorganFP_byId(self, id):
        return self.cur.execute('''SELECT morganfp FROM compounds 
                where id = ?''', (id,)).fetchone()

    def show_all_id_name_smile(self):
        return self.cur.execute('''SELECT id,name,smiles FROM compounds''').fetchall()

    def show_best_mtt_result(self, id):
        return self.cur.execute('''SELECT EC50 FROM compounds 
                        LEFT JOIN MTT ON compounds.id = MTT.compound_id
                        where compounds.id = ?
                        ORDER by EC50 desc

                        ''', (id,)).fetchone()


if __name__ == '__main__':
    db = DataBase('../external/resources/data.db')
    # db.add_solutor('Вода')
    # print(db.show_all_id_name_smile())
    # print(similiaryty_list_return('C1(=CC=C(C=C1)C=2N=NN(C2C#N)C3=C4N(C5=C3C(=C(C(=C5F)F)F)F)C=CC=C4)OC', db.show_all_id_name_smile()))
    # print(mol_mass_from_smiles('C1=CCCC2C(CC1)CCCCCC2'))
    print(db.show_best_mtt_result(3))