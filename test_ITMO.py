import sqlite3


connection = sqlite3.connect("test.db")
cursor = connection.cursor()

try:
    cursor.execute("""
    CREATE TABLE Car (
        Regist CHAR(10) PRIMARY KEY,
        CurN CHAR(6) NOT NULL,
        Region INTEGER NOT NULL,
        Brand VARCHAR(50) NOT NULL,
        Color VARCHAR(50) NOT NULL,
        Power INTEGER NOT NULL CHECK(Power>50),
        CarYear INTEGER NOT NULL,
        Mileage INTEGER,
        UNIQUE (CurN,Region)
    )
    """)
except:
    print('Already exist')

prompts = ["""INSERT INTO Car VALUES('9844720488', 'E340BT', 77, 'Lada Granta','Red', 87, 2017, 35)""",
           """INSERT INTO Car VALUES('6239572784', 'H109OK', 178, 'Volkswagen Polo','Blue', 105)""",
           """INSERT INTO Car VALUES('4752909757', 'A822EY', 99, 'Skoda Rapid','Black', 125, 2021, 35)""",
           """INSERT INTO Car VALUES('7984672834', 'T120AA', 98,'Black','Hyundai Solaris', 123, 2019, 20)""",
           """INSERT INTO Car VALUES('7478679847', 'B971HP', 199, 'Kia Sportage','White', 184, 2017, 35)""",
           """INSERT INTO Car VALUES('4728472878', 'T120AA', 77, 'Toyota RAV4','Silver-grey', 146, 2008)""",
           """INSERT INTO Car VALUES('4782487387', 'H454EE', 98, 'Skoda Rapid','Black', 75, 2021, 0)""",
           """INSERT INTO Car VALUES('4782487387', 'O638OA', 173, 'Mitsubishi Outlander','White', 230, 2021)""",
           """INSERT INTO Car VALUES('7284728297', 'M118EA', 66, 'Hyundai Solaris','Blue', 123, 2021, 0)""",
           """INSERT INTO Car VALUES('8779854025', 'A822EY', 99, 'Lada Granta','White', 87, 2017, 54)"""
           ]

# for i in range(0, 10):
#     try:
#         cursor.execute(prompts[i])
#         print(f'{i + 1} --- True')
#     except:
#         print(f'{i + 1} --- False')

cursor.execute(prompts[9])



connection.commit()
connection.close()

