import sys
sys.path.append('/usr/local/lib/python3.9/site-packages/selenium-3.141.0-py3.9.egg')
from selenium import webdriver

from selenium.webdriver.common.keys import Keys


chrome_options = webdriver.ChromeOptions()
chrome_options.add_argument(' — headless')
chrome_options.add_argument(' — no-sandbox')
chrome_options.add_argument(' — disable-dev-shm-usage')



driver = webdriver.Chrome('/Users/hongjianyang/chromedriver')
path = 'https://www.purpleair.com/sensorlist?key=ENF09HH3YHHVB8LW&show='
index = 31361
driver.get(path + str(index))

inp = driver.find_element_by_id("startdatepicker")
inp.send_keys("01/12/2021")

inp = driver.find_element_by_id("enddatepicker")
inp.send_keys("01/13/2021")

inp = driver.find_element_by_xpath("//*[@id='resultTable']/tbody/tr[1]/td[1]/input")
inp.click()

inp = driver.find_element_by_id("773512_download_button")
inp.send_keys(Keys.ENTER)

