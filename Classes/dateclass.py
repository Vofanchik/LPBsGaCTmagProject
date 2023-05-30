from datetime import datetime


def from_dot_to_rec(how, dt):
    if how == 0:
        return datetime.strptime(dt, '%d.%m.%Y').strftime('%Y-%m-%d')
    else:return datetime.strptime(dt, '%Y-%m-%d').strftime('%d.%m.%Y')

def from_str_to_date(how, date_str):
    if how == 0:
        return datetime.strptime(date_str, '%d.%m.%Y')
    else: return datetime.strptime(date_str, '%Y-%m-%d')