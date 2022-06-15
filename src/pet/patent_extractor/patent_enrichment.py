# -*- coding: utf-8 -*-

"""Script for extracting patent literature from SureChEMBL."""

import logging
import os
import time
from typing import Tuple

import pandas as pd
from tqdm import tqdm

from pemt.constants import DATA_DIR, PATENT_DIR, VALID_CODES

# Selenium specific settings
try:
    from selenium import webdriver
except ImportError:
    raise ValueError("please install selenium before running this script")

from selenium.webdriver.chrome.options import Options
from selenium.common.exceptions import NoSuchElementException

chrome_options = Options()
chrome_options.add_argument("--headless")
chrome_options.add_argument("--window-size=1920x1080")
chrome_options.add_argument("--disable-notifications")
chrome_options.add_argument('--no-sandbox')
chrome_options.add_argument('--verbose')
chrome_options.add_experimental_option("prefs", {
    "download.prompt_for_download": False,
    "download.directory_upgrade": True,
    "safebrowsing_for_trusted_sources_enabled": False,
    "safebrowsing.enabled": False
})
chrome_options.add_argument('--disable-gpu')
chrome_options.add_argument('--disable-software-rasterizer')

logger = logging.getLogger('__name__')
logger.setLevel(logging.INFO)

"""Constant factors related to scraping"""

os.makedirs(f'{PATENT_DIR}', exist_ok=True)


def get_valid_patent_list(
    schembl_id: str,
    system: str,
    chrome_driver_path: str,
    year: int
) -> Tuple[set, int]:
    """Get valid patents from SureChEMBL based on their IPC criteria and time period.

    :param schembl_id: The SureChEMBL id of the compound.
    :param system: The OS on which the code is running. It can be either of these: linux, mac, window.
    :param chrome_driver_path: The path of the chrome driver is located.
    :param year: The cutt-off year for searching the patent documents
    """

    # Replace path to chrome driver (https://sites.google.com/a/chromium.org/chromedriver/home)
    driver = webdriver.Chrome(
        options=chrome_options,
        executable_path=chrome_driver_path,
    )

    # function to take care of downloading file
    driver.command_executor._commands["send_command"] = ("POST", '/session/$sessionId/chromium/send_command')
    params = {'cmd': 'Page.setDownloadBehavior', 'params': {'behavior': 'allow', 'downloadPath': DATA_DIR}}
    driver.execute("send_command", params)

    system = system.lower()

    logger.debug('Getting page')
    driver.get(f'https://www.surechembl.org/chemical/{schembl_id}')
    logger.debug('Page done')

    time.sleep(8)

    # Go to patent tab
    try:
        patent_button = driver.find_element_by_xpath(
            '/html/body/div/div/div[2]/div/div/div[3]/div[2]/ul/li[3]'
        )
        patent_button.click()
    except NoSuchElementException:  # no patents found
        return set(), 0

    time.sleep(15)

    # Get the link for opening patent table
    try:
        new_link = driver.find_element_by_xpath(
            '/html/body/div/div/div[2]/div/div/div[3]/div[2]/div[3]/div[3]/a'
        ).get_attribute('href')
    except NoSuchElementException:  # no patents found
        return set(), 0

    # Get total number of patents
    range_val = int(driver.find_element_by_xpath(
        "//span[@class='total_hits_data']"
    ).text.replace(',', ''))

    driver.get(new_link)
    time.sleep(2)

    patent_info = set()

    logger.debug(f'Looking into {range_val} patents')

    for patent_count in range(1, range_val + 1):
        for i in range(2, 52):  # max number of elements in each page
            try:
                if system not in ['linux', 'mac']:
                    ipc_num = driver.find_element_by_xpath(
                        f'/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div[2]/table/tbody/tr[{i}]/td[4]/div[1]/table/tbody/tr/td[1]'
                    ).text
                else:
                    ipc_num = driver.find_element_by_xpath(
                        f'/html/body/div/div/div[2]/div[1]/div[2]/div/div[2]/table/tbody/tr[{i}]/td[4]/div[1]/table/tbody/tr/td[1]'
                    ).text
            except NoSuchElementException:  # No IPC code found
                continue

            if not ipc_num:  # cases where IPC number is missing are skipped
                continue

            # Filter based on codes
            code = ipc_num.split()[0]
            if code not in VALID_CODES:
                continue

            if system not in ['linux', 'mac']:
                patent_date = driver.find_element_by_xpath(
                    f'/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div[2]/table/tbody/tr[{i}]/td[3]'
                ).text
            else:
                patent_date = driver.find_element_by_xpath(
                    f'/html/body/div/div/div[2]/div[1]/div[2]/div/div[2]/table/tbody/tr[{i}]/td[3]'
                ).text
            patent_year = int(patent_date.split('-')[0])

            # Filter based on patent year
            if patent_year < year:
                continue

            if system not in ['linux', 'mac']:
                patent_number = driver.find_element_by_xpath(
                    f'/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div[2]/table/tbody/tr[{i}]/td[2]'
                ).text.split('\n')[1]
            else:
                patent_number = driver.find_element_by_xpath(
                    f'/html/body/div/div/div[2]/div[1]/div[2]/div/div[2]/table/tbody/tr[{i}]/td[2]'
                ).text.split('\n')[1]

            # Get assignee information
            if system not in ['linux', 'mac']:
                assignee = driver.find_element_by_xpath(
                    f'/html/body/div[1]/div/div[2]/div[1]/div[2]/div/div[2]/table/tbody/tr[{i}]/td[4]/div[1]/table/tbody/tr/td[2]/a'
                ).text
            else:
                assignee = driver.find_element_by_xpath(
                    f'/html/body/div/div/div[2]/div[1]/div[2]/div/div[2]/table/tbody/tr[{i}]/td[4]/div[1]/table/tbody/tr/td[2]/a'
                ).text
            patent_info.add((patent_number, patent_date, ipc_num, assignee))

        if range_val < 50:  # entry fits in 1 page
            break

        # Go to next page
        try:
            if patent_count == 1:
                nx_button_num = 2
            else:
                nx_button_num = 4
            if system not in ['linux', 'mac']:
                next_page = driver.find_element_by_xpath(
                    f'/html/body/div[1]/div/div[2]/div[1]/div[2]/div[1]/div[3]/div[2]/ul/li[{nx_button_num}]/a'
                ).get_attribute('href')
            else:
                next_page = driver.find_element_by_xpath(
                    f'/html/body/div/div/div[2]/div[1]/div[2]/div[1]/div[3]/div[2]/ul/li[{nx_button_num}]/a'
                ).get_attribute('href')
            driver.get(next_page)
            time.sleep(8)
        except NoSuchElementException:
            continue

    return patent_info, range_val


def extract_patent(
    analysis_name: str,
    chrome_driver_path: str,
    os_system: str = 'linux',
    patent_year: int = 2000
) -> pd.DataFrame:
    """Extract and store all valid patent document metadata.

    :param analysis_name: Name of the analysis.
    :param os_system: The OS on which the code is running. It can be either of these: linux, mac, window.
    :param chrome_driver_path: The path of the chrome driver is located.
    """
    df = pd.read_csv(
        f'{PATENT_DIR}/{analysis_name}_chemicals.tsv',
        sep='\t',
        dtype=str,
        usecols=[
            'chembl',
            'schembl_id'
        ]
    )

    if df.empty:
        return pd.DataFrame()

    os_system = os_system.lower()
    assert os_system in ['linux', 'mac', 'windows']
    logger.warning(f'Currently running on {os_system} OS. Please change if this is not the case.')

    # Check for existing cache file
    if os.path.exists(f'{PATENT_DIR}/{analysis_name}_patent_data.tsv'):
        patent_df = pd.read_csv(f'{PATENT_DIR}/{analysis_name}_patent_data.tsv', sep='\t')
    else:
        patent_df = pd.DataFrame(columns=['chembl', 'surechembl'])

    cache_count = 0

    for chembl_id, surechembl_idx in tqdm(df.values, total=df.shape[0]):

        if pd.isna(surechembl_idx):
            continue

        _info_df = patent_df[
            (patent_df['chembl'] == chembl_id) & (patent_df['surechembl'] == surechembl_idx)
        ]

        if not _info_df.empty:
            continue

        cache_count += 1

        patent_info, total = get_valid_patent_list(
            schembl_id=surechembl_idx,
            system=os_system,
            chrome_driver_path=chrome_driver_path,
            year=patent_year
        )

        if len(patent_info) == 0:
            patent_df = pd.concat(
                [
                    patent_df,
                    pd.DataFrame({
                        'chembl': chembl_id,
                        'surechembl': surechembl_idx,
                        'patent_id': '',
                        'date': '',
                        'ipc': '',
                        'assignee': ''}, index=[0])
                ],
                ignore_index=True
            )
        else:
            for patent in patent_info:
                (
                    pid, date, ipc_code, assignee
                ) = patent
                patent_df = pd.concat(
                    [
                        patent_df,
                        pd.DataFrame({
                            'chembl': chembl_id,
                            'surechembl': surechembl_idx,
                            'patent_id': pid,
                            'date': date,
                            'ipc': ipc_code,
                            'assignee': assignee}, index=[0])
                    ],
                    ignore_index=True
                )

        if cache_count == 5:  # in case of internet issues
            patent_df.drop_duplicates(inplace=True)
            patent_df.to_csv(f'{PATENT_DIR}/{analysis_name}_patent_data.tsv', sep='\t', index=False)
            cache_count = 0

    patent_df.drop_duplicates(inplace=True)
    patent_df.to_csv(f'{PATENT_DIR}/{analysis_name}_patent_data.tsv', sep='\t', index=False)
    return patent_df
